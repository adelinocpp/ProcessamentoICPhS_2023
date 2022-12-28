#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 02:18:32 2022

@author: adelino
"""
import os
from pathlib import Path
import numpy as np
from scipy.fft import fft
from scipy.integrate import simpson
from chardet.universaldetector import UniversalDetector
import subprocess
# -----------------------------------------------------------------------------
def simpson_integral(t,f):
    Nf = len(f)
    N = len(t)
    if not (N == Nf):
        print('Inetgral error: diffent number of points')
        return 0
    h = (t[-1]-t[0])/(N - 1)
    I_simp = (h/3) * (f[0] + 2*np.sum(f[:N-2:2]) \
            + 4*np.sum(f[1:N-1:2]) + f[N-1])
    return I_simp
# -----------------------------------------------------------------------------
def spectral_ratios(audio, sr, time_step, nFFT = 1024):
    LTF = []
    COG_1 = []
    COG_2 = []
    COG_23 = []
    nStep = int(time_step*sr)
    if (nStep > nFFT):
        nFFT = int(2**np.ceil(np.log2(nStep)))
    hFFT = int(0.5*nFFT)
    f = np.linspace(0,0.5*sr,hFFT)
    wAudio = np.zeros((nFFT,))
    lowFilter = np.zeros((len(f)))
    idxLow = (f <= 850).nonzero()[0]
    idxSmt = np.multiply((f >= 850),(f <= 950)).nonzero()[0]
    lowFilter[idxLow] = 1
    lowFilter[idxSmt] = (950 - f[idxSmt])/100
    for i in range(0,len(audio)-nStep,nStep):
        
        wAudio[:nStep] = audio[i:(i+nStep)]
        faudio = np.abs(fft(wAudio,n=nFFT)[:hFFT])
        lfaudio = np.multiply(lowFilter,faudio)
        lfpower = simpson(lfaudio,x=f)
        
        fpower_1 = simpson(faudio,x=f)
        fpower_2 = simpson(np.power(faudio,2),x=f)
        fpower_23 = simpson(np.power(faudio,2/3),x=f)
        fint_1 = simpson(np.multiply(f,faudio),x=f)
        fint_2 = simpson(np.multiply(f,np.power(faudio,2)),x=f)
        fint_23 = simpson(np.multiply(f,np.power(faudio,2/3)),x=f)
        COG_1.append(fint_1/fpower_1)
        COG_2.append(fint_2/fpower_2)
        COG_23.append(fint_23/fpower_23)
        LTF.append(2*(fpower_1 - lfpower)/sr)
    return COG_1, COG_2, COG_23, LTF
# -----------------------------------------------------------------------------
def list_contend(folder='./', pattern=()):
    pattern = tuple([x.upper() for x in pattern])
    list_files = [os.path.join(root, name)
             for root, dirs, files in os.walk(folder)
                 for name in files
                     if name.upper().endswith(pattern)]
    list_files.sort()
    return list_files
# -----------------------------------------------------------------------------
def convert_utf8(filename,codefrom):

    tempFilename = '{:}/tempfilename.TextGrid'.format(Path(filename).parent)
    # renameComand = "mv \"{:}\" \"{:}\"".format(os.path.abspath(filename),os.path.abspath(tempFilename))
    convertcomand = 'iconv -f {:} -t utf-8 {:} -o {:}'.format(codefrom,os.path.abspath(tempFilename),os.path.abspath(filename))
    deleteComand = 'rm {:}'.format(os.path.abspath(tempFilename))
    if os.path.exists(filename):
        os.rename(os.path.abspath(filename),os.path.abspath(tempFilename))
        # TODO: problema na conversao
#        returnStr =  subprocess.Popen(cmdString, shell=True, stdout=subprocess.PIPE).stdout
#        self.codec =  returnStr.read().decode()[0:-1]
        subprocess.Popen(convertcomand, shell=True, stdout=subprocess.PIPE).stdout
        subprocess.Popen(deleteComand, shell=True, stdout=subprocess.PIPE).stdout
        
# -----------------------------------------------------------------------------
def check_utf8(filename):
    resutl_utf8 = False
    oriCode = ''
    if os.path.exists(filename):
        u = UniversalDetector()
        u.reset()
        with open(filename, "rb") as f:
            for index, line in enumerate(f):
                u.feed(line)
                if u.done:
                    break
        u.close()
        oriCode = u.result['encoding'].lower()
        resutl_utf8 = oriCode == 'utf-8'
    return resutl_utf8, oriCode
# -----------------------------------------------------------------------------
def textgrid_to_interval_matrix(filename, anyLabel=1, labelValue="", convertNumber=0,tierNumber=-1):
    '''filename:   endereço completo do arquivo .textgrid
       anyValue:   Se anyLabel = 0 lista os intervalos indicados em labelValue
                   Se anyLabel = 1 (padrão) lista apenas os intervalos com etiqueta diferente de vazio
                   Se anyLabel = -1 lista todos intervalos inclusive os vazios
       labelValue: valor da etiqueta (do inervalo) a ser buscado. Padrão = "". Se especificado
                   retorna apenas os intervalos que contem a etiqueta igual a labelValue
       convertNumber: Se convertNumber=0 (padrão) retorna o valor da etiqueta como string
                      Se convertNumber=1 tenta realizar a conversão do conteúdo da etiqueta para valor numérico
       tierNumber:    Especifica a camada para realizar a busca.
                     Se tierNumber=-1 (padrão) realiza a busca em todas as camadas
                     A primeira camada é tierNumber=0
    '''
    tierMatrix = []
    file_suffix = Path(filename).suffix.upper()
    if (file_suffix != ".TEXTGRID"):
        print("Erro: arquivo não possui extensão do tipo textgrid.")
    
    # Verifica e converte para UTF8
    
    
    try:
        isUTF8, fromCode = check_utf8(filename)
        if (not isUTF8):
            print('Nao e UFT-8')
            convert_utf8(filename,fromCode)
        with open(filename, 'r') as file:
            # --- Faz a leitura do arquivo como uma lista de linhas
            TextGridLines = file.readlines()
        
        fileTypeStr = TextGridLines[0].split("= ")[1]
        # objClassStr = TextGridLines[1].split("= ")[1]
        # xmin = float(TextGridLines[3].split("= ")[1])
        # xmax = float(TextGridLines[4].split("= ")[1])
        tiersBool = TextGridLines[5].split("? ")[1]
        nTiers = int(TextGridLines[6].split("= ")[1])
        
        if (fileTypeStr != '"ooTextFile"\n'):
            print("Erro: Arquivo textgrid indica que não é do tipo TextFile.")
            return tierMatrix
        if (fileTypeStr != '"ooTextFile"\n'):
            print("Erro: Arquivo textgrid indica conteúdo não é do tipo TextGrid.")
            return tierMatrix
        if (tiersBool != '<exists> \n'):
            print("Erro: Arquivo textgrid indica que não possui camadas.")
            return tierMatrix
    except:
        print("Erro:Problema de leitura do arquivo textgrid")  
        return tierMatrix
    
    if (tierNumber > 0) and (tierNumber >= nTiers):
        print("Erro: Indicadca camada {:}. O arquivo textgrit possui {:} camadas. Selecione uma camada entre 0 e {:}.".format(tierNumber,nTiers,nTiers-1))
        
    
    initTiers = []
    kTier = 1
    for idx, line in enumerate(TextGridLines):
        tierNumberStr = '    item [{:}]:\n'.format(kTier)
        if (line == tierNumberStr):
            initTiers.append(idx)
            kTier += 1
    
    tierMatrix = []
    for idxT, initLine in enumerate(initTiers):
        strClass = TextGridLines[initLine + 1].split("= ")[1]
        if (strClass != '"IntervalTier" \n'):
            print("Aviso: Saltando camada {:} que não é do tipo de intervalos...".format(idxT))
            
            nIntervals = int(TextGridLines[initLine + 5].split("= ")[1])
            newTier = []
            for i in range(initLine + 6,initLine + 6 + 3*nIntervals,3):
                itv_xmin = float(TextGridLines[i+1].split("= ")[1])
                itv_textT = TextGridLines[i+2].split("= ")[1].strip().replace('"','')
                insertValue = True
                if (anyLabel == 0) and (itv_textT != labelValue): 
                    insertValue = False  
                if (anyLabel == 1) and (len(itv_textT) == 0):
                    insertValue = False     
                if insertValue:
                    if (convertNumber):
                        try:
                            itv_text = float(itv_textT)
                        except:
                            print("Erro: Conversao para valor numerico não e possivel de valor {:}".format(itv_textT))
                    else:
                        itv_text = itv_textT
                    newTier.append([itv_xmin,itv_text])  
        else:
            nIntervals = int(TextGridLines[initLine + 5].split("= ")[1])
            newTier = []
            for i in range(initLine + 6,initLine + 6 + 4*nIntervals,4):
                itv_xmin = float(TextGridLines[i+1].split("= ")[1])
                itv_xmax = float(TextGridLines[i+2].split("= ")[1])
                itv_textT = TextGridLines[i+3].split("= ")[1].strip().replace('"','')
                insertValue = True
                if (anyLabel == 0) and (itv_textT != labelValue): 
                    insertValue = False  
                if (anyLabel == 1) and (len(itv_textT) == 0):
                    insertValue = False     
                if insertValue:
                    if (convertNumber):
                        try:
                            itv_text = float(itv_textT)
                        except:
                            print("Erro: Conversao para valor numerico não e possivel de valor {:}".format(itv_textT))
                    else:
                        itv_text = itv_textT
                    newTier.append([itv_xmin,itv_xmax,itv_text])    
        tierMatrix.append(newTier)    
        
    return tierMatrix  
            
            
            
            
            
            
            
            
            
            
            

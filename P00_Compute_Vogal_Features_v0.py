#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cálculo de características acústicas em vogais
Minha humilde contribuição para o trabalho de Maria Cantoni e Thais Cristofaro

---
Data: 09/11/2022
Autor: Adelino P. Silva
mail: adelinocpp@gmail.com, adelinocpp@gmail.com, adelino@ufmg.br
tel: +55 31 98801-3605

--- 
Rotina para a leitura de arquivos do tipo .TextgGrid (praat) e calculo de 
taxas espectrais como:
1 - HNR (de Signal_Analysis.features.signal)
2 - centro de gravidade espectral do tipo 1, 2 e 2/3, (vide link abaixo) e

Centro de gravidade praat:
https://www.fon.hum.uva.nl/praat/manual/Spectrum__Get_centre_of_gravity___.html

----
Contribuições basicas de funções:
textgrid_to_interval_matrix: decompõe o arquivo .Textgrid em uma matriz de intervalos
spectral_ratios: calcula os centros de gravidade e LTF
---
Livre para uso e modificações.
Em caso de dúvidas entre em contato.
"""
from utils.file_utils import list_contend, textgrid_to_interval_matrix, spectral_ratios
from scipy.io import wavfile
from Signal_Analysis.features.signal import get_HNR
import numpy as np
from pathlib import Path
import sys
import os
# import praat_formants_python as pfp
from utils.formant_lpc import format_lpc, intensity
from g2p.g2p import G2PTranscriber
import re
import warnings
warnings.filterwarnings("ignore")
# -----------------------------------------------------------------------------
def getMeanPercentualInterval(variable,pIni,pFim):
    nPoints = len(variable)
    nIni = int(np.floor(pIni*nPoints))
    nFim = int(np.ceil(pFim*nPoints))
    return np.mean(variable[nIni:nFim])
# -----------------------------------------------------------------------------
def is_ditongo(strVogal):
    # vogal = ('a','e','i','o','u','A','E','I','O','U')
    vogal_ext = ('a','e','i','o','u','A','E','I','O','U','ã','â','à','á',
                 'Ã','Â','À','Á','é','ê','É','Ê','í','Í','ó','õ','ô','Ó','Õ','Ô',
                 'ú','ü','Ú','Ü')
    strRes = ''
    for i in strVogal:
        if (i in vogal_ext):
            strRes = strRes + i
    return int(len(strRes) > 1)
# -----------------------------------------------------------------------------
def tag_to_bass_vowel(tag, idTag):
    listConsoante = ['b','c','d','f','g','h','j','k','l','m','n','p','q','r',
                     's','t','v','w','x','y','z']
    vogal_bas = np.array(['a','e','i','o','u'])
    # vogal_gra = np.array(['a','ã','â','à','á','e','é','ê','i','í','ĩ','o','ó','õ','ô','u','ú','ü','ũ'])
    # vogal_idx = np.array([0,0,0,0,0,1,1,1,2,2,2,3,3,3,4,4,4,4])
    vogal_gra = np.array(["a","á","à","â","ɐ","ã","e","é","ɛ","ê","i","ɪ","ɪ̃","í","ĩ","o","ó","ɔ","ô","õ","u","ú","ü","ʊ","ũ"])
    vogal_idx = np.array([0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4])
    newTag = ''
    for lk in tag:
        if not (lk in listConsoante):
            newTag = newTag + lk
    newTag = newTag.lower()
    try:
        idxg = (newTag == vogal_gra).nonzero()[0][0]
        idxb = vogal_idx[idxg]
        return vogal_bas[idxb], idxb
    except:
        print("Erro convert vowel: IDX {:}, {:}. ID: {:}".format(idxg,idxb,idTag))
        return '0', -1
# -----------------------------------------------------------------------------
# percorre uma silaba e tenta identificar se existe uma vogal
def has_vogal(phono):
    vogal_gra = np.array(["ɐ","a","á","à","â","ɐ","ã","e","é","ɛ","ê","i","ɪ","ɪ̃","í","ĩ","ĩ","o","ó","ɔ","ô","õ","u","ú","ü","ʊ","ũ"])
    ret_val = False
    for g in phono:
         ret_val = ret_val or (g in vogal_gra)
    return ret_val
# -----------------------------------------------------------------------------
def pos_indicated_vowel(phono, vowel, idTag):
    # vogal_bas = np.array(['a','e','i','o','u'])
    
    # vogal_gra = ('a','ã','â','à','á','e','é','ê','i','í','o','ó','õ','ô','u','ú','ü')
    vogal_pho = np.array(["a","á","à","â","ɐ","ã","e","é","ɛ","ê","i","ɪ","í","ĩ","ɪ̃","o","ó","ɔ","ô","õ","u","ú","ü","ʊ","ũ"])
    tabBasPhono = np.array([[0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1, 1, 1, 1, 3 ],
                            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 0, 0, 0, 0, 3 ]])
        
    v_pos = np.array([])
    v_dist = np.array([])
    
    # for i, g in enumerate(graph):
    #     if g in vogal_pho:
    #         g_bas, i_bas = tag_to_bass_vowel(g,idTag)
    #         try:
    #             i_pho = (g == vogal_pho).nonzero()[0][0]
    #             v_pos = np.append(v_pos, i)
    #             v_dist = np.append(v_dist, tabBasPhono[i_bas,i_pho] )
    #         except:
    #             print("Erro find vowel {:}. ID: {:}".format(g,idTag))
            
    g_bas, i_bas = tag_to_bass_vowel(vowel,idTag)
    for i, g in enumerate(phono):
        if g in vogal_pho:
            try:
                i_pho = (g == vogal_pho).nonzero()[0][0]
                v_pos = np.append(v_pos, i)
                v_dist = np.append(v_dist, tabBasPhono[i_bas,i_pho])
            except:
                print("Erro find vowel {:}. ID: {:}".format(g,idTag))
            
    u_pos = np.unique(v_pos)
    if (len(u_pos) == 1):
        return int(u_pos[0])
    elif (len(u_pos) == 0):
        print("Erro: ID: {:}".format(idTag))
    else:
        u_dist = 5*np.ones(u_pos.shape)
        for i, u in enumerate(u_pos):
            idx = (u == u_pos).nonzero()[0]
            u_dist[i] = np.mean(v_dist[idx])
        idx = np.argmin(u_dist)
        return int(u_pos[idx])
# -----------------------------------------------------------------------------
def estimate_syllabe_position(phonografico, numSyl, posVogal):
    returnValue = -1
    if (len(phonografico) != numSyl):
        return returnValue
    
    # posicao da tónica do final da palavra (direita) para i início (esquerda)
    stressPos = numSyl - np.array(["ˈ" in parte[j] for parte in phonPalavra]).nonzero()[0][0] - 1
    
    if(stressPos < 3):
        if (posVogal == 0):
            returnValue = (numSyl - stressPos - 1)
        # Explicação acima
        elif (posVogal < 3):
            returnValue = (numSyl - 1 - (2-posVogal))
        elif (posVogal > 2):    
            returnValue = (numSyl - posVogal + 1 - stressPos)    
    else:
        print("9: Sílaba tónica em posição {:}, diferente da esperada.".format(stressPos))
    
    return returnValue
    '''
    if (stressPos == 0):
        if (posVogal == 0):
            returnValue = (numSyl - stressPos - 1)
        
    elif(stressPos == 1):
        if (posVogal == 0):
            returnValue = (numSyl - stressPos - 1)
        # Para pos tonica:
        # Exemplos   len   posVogal 2(ultimo no vetor)   posVogal 1(penultimo no vetor)
        #    3012:    4    4-1 = 3 [len-1 -(2-posVogal)] 4-2 = 2 [len-1 -(2-posVogal)]
        #     302:    3    3-1 = 2 [len-1 -(2-posVogal)]  X: não tem, não ocorre
        elif (posVogal < 3):
            returnValue = (numSyl - 1 - (2-posVogal))
        # Para pre tonica:
        elif (posVogal > 2):    
            returnValue = (numSyl - posVogal + 1 - stressPos)
            
    elif(stressPos == 2):
        if (posVogal == 0):
            returnValue = (numSyl - stressPos - 1)
        # Explicação acima
        elif (posVogal < 3):
            returnValue = (numSyl - 1 - (2-posVogal))
        elif (posVogal > 2):    
            returnValue = (numSyl - posVogal + 1 - stressPos)    
    else:
        print("9: Sílaba tónica em posição {:}, diferente da esperada.".format(stressPos))
    
    '''
# -----------------------------------------------------------------------------


AUDIO_FOLDER = '../Audios/'
CSVFILE = './csvDataAudios.csv'
audiofiles = list_contend(folder=AUDIO_FOLDER, pattern=('.wav',))
textgridfiles = list_contend(folder=AUDIO_FOLDER, pattern=('.textgrid',))

if (len(audiofiles) != len(textgridfiles)):
    print("Erro: número de arquivos de áudio não corresponde ao numero de TextGrid")

# Tamanho do passo de tempo em segundos
valStep = 0.005
valWin  = 0.020
listConsoante = ('b','c','d','f','g','h','j','k','l','m','n','p','q','r','s','t','v','w','x','y','z')

list_phon_vowels = ("a", "e", "o", "á", "é", "í", "ó", "ú", "ã", "õ", "â", "ê", "ô", "à", "ü","ɪ", "ɛ", "ɐ", "ʊ", "ɔ")

list_phon_consonant = ('b','c','d','f','g','h','j','k','l','m','n','p','q','r','s','t','v','w','x','y','z','ʃ', 'ʎ', 'ɲ', 'ɳ', 'ɾ', 'ɣ', 'ʒ', 'ʤ', 'ʧ', 'ŋ')

csvLines = []
strTitle = "ID, Duração, F1, F2, F1_b, F2_b, intensidade, tag, Posiçao, HNR, Fonetica, Fonologica, Ditongo, Palavra, Tonicidade, Precedente, Seguinte, Fechada\n".replace(",","\t")
csvLines.append(strTitle)
useTiers = (0,)
tabId = 0;
maxNSyllab = 0

for idxtg, tgFile in enumerate(textgridfiles):
    value = textgrid_to_interval_matrix(tgFile)    
    sr, audio = wavfile.read(audiofiles[idxtg])
    if (len(audio.shape) > 1):
        nSamples, nChannel = audio.shape
    else:
        nSamples  = len(audio)
        nChannel  = 1
    if (nChannel > 1):
        audio = np.mean(audio,axis = 1)
    audio = audio/np.max(np.abs(audio))
        
    for j in useTiers:
        intervalMtx = value[j]    
        for idxL, interval in enumerate(intervalMtx):
            if (tabId == 4):
                print('Depurando...')
            tabDitongo = 0
            nIni = int(interval[0]*sr)
            nFim = int(interval[1]*sr)
            tabDuration = interval[1] - interval[0]
            if (tabDuration < (valWin + valStep)):
                print("1: Segmento ID {:} de {:} seg é muito pequeno para processamento.".format(tabId,tabDuration))
                print("1: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
                
            tags = interval[2].split("-")
            
            if not (len(tags) == 5):
                print("2: Etiqueta {:} ID {:} com problema no número de marcações.".format(interval[2],tabId))
                print("2: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
            if (int(tags[2]) > maxNSyllab):
                maxNSyllab = int(tags[2])
                
            tabFonetica = tags[0]
            if (len(tabFonetica) == 1) and (tabFonetica.lower() in listConsoante):
                print("2: Etiqueta {:} ID {:} apresenta tag que nao e vogal {:}.".format(interval[2],tabId,tabFonetica))
                continue
            
            tabPalavra = tags[1]
            tabTonicidade = tags[3]
            tabDitongo = is_ditongo(tags[0])
            # if (len(tags[0]) == 1):
            #     tabDitongo = 0
            # elif(len(tags[0]) == 2):
            #     tabDitongo = 1
            # else:
            #     print('3: Tag {:} ID {:} nao se encaixa em vogal ou ditongo!'.format(tags[0],tabId))
            #     print("3: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                
            selAudio = audio[nIni:nFim]
            form2, _ = format_lpc(selAudio,sr,winstep=valStep,winlen=valWin)
            inten = intensity(selAudio,sr,winstep=valStep,winlen=valWin)
            tabMeanHNR, _ = get_HNR(selAudio,sr,time_step=valStep,periods_per_window = 1.875)
            tabIntensity = getMeanPercentualInterval(inten,0.2,0.8)
            if (tabDitongo == 0):
                tabF1 = getMeanPercentualInterval(form2[0,:],0.2,0.8)
                tabF2 = getMeanPercentualInterval(form2[1,:],0.2,0.8)
                tabF1_b = -1
                tabF2_b = -1
            else:
                tabF1 = getMeanPercentualInterval(form2[0,:],0.1,0.2)
                tabF2 = getMeanPercentualInterval(form2[1,:],0.1,0.2)
                tabF1_b = getMeanPercentualInterval(form2[0,:],0.8,0.9)
                tabF2_b = getMeanPercentualInterval(form2[1,:],0.8,0.9)
            
            # momento da transcricao groafica para fonetica
            g2p = G2PTranscriber(tabPalavra, algorithm='ceci')
            phonPalavra = g2p.transcriber().split(",")[0].split('.')
            grafPalavra = g2p.get_syllables_with_hyphen().split('-')
            
            if (len(grafPalavra) != len(phonPalavra)):
                phonPalavra = phonPalavra[:(len(grafPalavra)-1)]
            
            if (len(grafPalavra) != int(tags[2])):
                print('4: ID {:}. Problema na silabificação (tag {:}) de \"{:}\" como {:} , {:}'.format(tabId,int(tags[2]),tabPalavra,phonPalavra, grafPalavra))
                print("4: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
            
            if (int(tags[4]) >  len(grafPalavra)):
                print("5: Etiqueta {:} ID {:} com problema. Possicao maior que o numero de silabas".format(interval[2],tabId))
                print("5: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
            
            vogalPos = estimate_syllabe_position(phonPalavra,int(tags[2]),int(tags[4]))
            if (vogalPos < 0):
                print("6: Etiqueta {:} ID {:} com problema. Falha na estimaçao da posiçao da silaba".format(interval[2],tabId))
                print("6: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
            
            # remove o marcador de silaba tonica
            phonSilaba = phonPalavra[vogalPos].replace("ˈ","")
            grapSilaba = grafPalavra[vogalPos]
            tabFonologico = phonSilaba #phonPalavra[vogalPos]
            # TODO:
                # Nao identifica algumas vogais tipo em " ʧĩ"
                # nasais aparecem com tamanho 2, tipo "ɪ̃"
            
            hasNasal = False
            idxNasal = []
            if ('͂' in phonSilaba) or ('̃' in phonSilaba):
                hasNasal = True    
                oriphonSilaba = phonSilaba
                # nasalIdx = phonSilaba.find('͂')
                phonSilaba = phonSilaba.replace('͂','').replace('̃','')
                print("Tratar as nasais")
                
            
            if (not has_vogal(tabFonologico)):
                print("7: ID {:} com problema, {:} sem vogal.".format(tabId, tabFonologico))
                continue
            
            tabPrecedente = 'NA'
            tabSeguinte = 'NA'
            if (tabDitongo):
                vogLen = len(grapSilaba)-1
                pos = grapSilaba.find(tags[0])
                if (pos == 0) and (len(phonSilaba) == vogLen+1):
                    tabPrecedente = 'NA'
                    tabSeguinte = 'NA'
                if (pos == 0) and (len(phonSilaba) > vogLen+1):
                    tabPrecedente = 'NA'
                    tabSeguinte = phonSilaba[1]
                if (pos > 0) and (len(phonSilaba) == (pos+vogLen+1)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = 'NA'
                if (pos > 0) and (len(phonSilaba) > (pos+vogLen+1)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = phonSilaba[pos+1]
            else:
                pos = pos_indicated_vowel(phonSilaba, tags[0], tabId)
                # pos = grapSilaba.find(tags[0])
                if (pos == 0) and (len(phonSilaba) == 1):
                    tabPrecedente = 'NA'
                    tabSeguinte = 'NA'
                if (pos == 0) and (len(phonSilaba) > 1):
                    tabPrecedente = 'NA'
                    tabSeguinte = phonSilaba[1]
                if (pos > 0) and (len(phonSilaba) == (pos+1)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = 'NA'
                if (pos > 0) and (len(phonSilaba) > (pos+1)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = phonSilaba[pos+1]
            
            if (hasNasal):
                print("recolocar a nasal")
                phonSilaba = oriphonSilaba
                
            tabFechada = int((grapSilaba[-1] in listConsoante))
            tabData = (tabId,tabDuration, tabF1,tabF2,tabF1_b,tabF2_b,tabIntensity,
                       interval[2],vogalPos,
                       tabMeanHNR, tabFonetica,tabFonologico, tabDitongo, 
                       tabPalavra, tabTonicidade, tabPrecedente, tabSeguinte,tabFechada)
            
            strData = "{:}\n".format(tabData).replace("(","").replace(")","").replace(",","\t")              
            # COG_1, COG_2, COG_23, LTF = spectral_ratios(selAudio,sr,time_step=valStep)
            # sys.exit("Saida de depuraçao")
            csvLines.append(strData)
            tabId = tabId + 1
        # sys.exit("Saida de depuraçao - RODOU APENAS CAMADA 1!")
        
with open(CSVFILE, 'w') as file:
        file.writelines(csvLines)

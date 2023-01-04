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
from unidecode import unidecode
warnings.filterwarnings("ignore")

abdist = np.genfromtxt ('dst_alfabeto.csv', delimiter=",")
ab_list=[]
for i in range(97,123):
    ab_list.append(chr(i))
# -----------------------------------------------------------------------------
'''
    Retorna a posiçao da vogal marcada na palavra
''' 
def pos_vowel_in_word(grafPalavra, vogalPos, mTag):
    preLetters = 0
    for i in range(0,vogalPos):
        preLetters = preLetters + len(grafPalavra[i])
    pos, valT = find_pos_of_tag(grafPalavra[vogalPos],mTag)
    return preLetters + pos
# -----------------------------------------------------------------------------
''' 
Calcula a "distancia" (diferença) entre letras do alfabeto em palavra
com base na tabela do arquivo "dst_alfabeto.csv".
ajusta algumas semelhanças como:
    dist_of_letters('u','u') = 0
    dist_of_letters('u','w') = 0.5
    dist_of_letters('u','l') = 0.5
    dist_of_letters('w','l') = 0.75
    dist_of_letters('i','e') = 0.5
    dist_of_letters('u','a') = 1
'''
def dist_of_letters(L1,L2):
    idxL1 = ab_list.index(unidecode(L1).lower())
    idxL2 = ab_list.index(unidecode(L2).lower())
    return abdist[idxL1,idxL2]
# -----------------------------------------------------------------------------
'''
Calcula a media de uma serie no tamanho percentual do intervalo. 
Se uma serie x possui tem 50 amostras e deseja-se calcular a media entre o 
ponto 5 e o 40 utiliza-se:
    getMeanPercentualInterval(x,0.1,0.8)
'''
def getMeanPercentualInterval(variable,pIni,pFim):
    nPoints = len(variable)
    nIni = int(np.floor(pIni*nPoints))
    nFim = int(np.ceil(pFim*nPoints))
    return np.mean(variable[nIni:nFim])
# -----------------------------------------------------------------------------
'''
Recebe uma marcacao que pode ser uma unica vogal ou um ditongo e retorna 
um valor logico indicando se a marcaçao e ditongo ou nao.
'''
def is_ditongo(strVogal):
    vogal_ext = ('a','e','i','o','u','A','E','I','O','U','ã','â','à','á',
                 'Ã','Â','À','Á','é','ê','É','Ê','í','Í','ó','õ','ô','Ó','Õ','Ô',
                 'ú','ü','Ú','Ü')
    strRes = ''
    for i in strVogal:
        if (i in vogal_ext):
            strRes = strRes + i
    return int(len(strRes) > 1)
# -----------------------------------------------------------------------------
'''
    Recebe como entrada um caracter de vogal com acentos e retorna o 
    caracter sem acento e index de identificaçao.
'''
def tag_to_bass_vowel(tag, idTag):
    listConsoante = ['b','c','d','f','g','h','j','k','l','m','n','p','q','r',
                     's','t','v','w','x','y','z']
    vogal_bas = np.array(['a','e','i','o','u'])
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
        print("Erro convert vowel: TAG {:}. ID: {:}".format(tag,idTag))
        return '0', -1
# -----------------------------------------------------------------------------
'''
    Encontra a posiçao de uma marcaçao em um trecho ortografico com base na 
    "distancia" entre as letras.
'''
def find_pos_of_tag(mSilabe,mTag):
    mSilabe = unidecode(mSilabe).lower()
    mTag = unidecode(mTag).lower()
    nSil = len(mSilabe)
    nTag = len(mTag)
    if (nTag > nSil):
        print("Problema com tamanho da marcaçao Silaba: {:}, TAG {:}".format(mSilabe,mTag))
        return -1, 1
    vDist = np.zeros((nSil-nTag+1,))
    for i in range(0,nSil-nTag+1):
        iDist = 0
        for j in range(0,nTag):
            iDist = iDist + dist_of_letters(mSilabe[i+j],mTag[j])/nTag
        vDist[i] = iDist
    posT = np.argmin(vDist)
    valT = np.min(vDist)
    return posT, valT
# -----------------------------------------------------------------------------
'''
    percorre uma silaba em transcriçao fonetica e tenta identificar se existe uma vogal
'''
def has_vogal(phono):
    vogal_gra = np.array(["ɐ","a","á","à","â","ɐ","ã","e","é","ɛ","ê","i","ɪ","ɪ̃","í","ĩ","ĩ","o","ó","ɔ","ô","õ","u","ú","ü","ʊ","ũ"])
    ret_val = False
    for g in phono:
         ret_val = ret_val or (g in vogal_gra)
    return ret_val
# -----------------------------------------------------------------------------
def pos_indicated_vowel(phono, vowel, idTag):
    vogal_pho = np.array(["a","á","à","â","ɐ","ã","e","é","ɛ","ê","i","ɪ","í","ĩ","ɪ̃","o","ó","ɔ","ô","õ","u","ú","ü","ʊ","ũ"])
    tabBasPhono = np.array([[0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ],
                            [1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1, 1, 1, 1, 3 ],
                            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 0, 0, 0, 0, 3 ]])
        
    v_pos = np.array([])
    v_dist = np.array([])
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
'''
Estima a posiçao da silaba na palavra, a partir do inicio (a esquerda) 
com base na seguinte codificaçao:
    0: tonica
    1: pos tonica
    2: pos pos tonica
    3: 1a pre-tonica
    4: 2a pre-tonica
    ....
'''
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
strTitle = "ID, Duração, F1, F2, F1_b, F2_b, intensidade, tag, Posiçao, HNR, Fonetica, Fonologica, Ditongo, Palavra, Tonicidade, Precedente, Seguinte, Fechada, Silabas, Oral, LetraPre, LetraSeg, Sexo, Arquivo\n".replace(",","\t")
csvLines.append(strTitle)
useTiers = (0,)
tabId = 0;
maxNSyllab = 0

for idxtg, tgFile in enumerate(textgridfiles):
    value = textgrid_to_interval_matrix(tgFile)    
    sr, audio = wavfile.read(audiofiles[idxtg])
    tabSexo = tgFile[-15]
    tabFileName = tgFile.split("/")[-1].split(".")[0]
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
            tabDitongo = 0
            nIni = int(interval[0]*sr)
            nFim = int(interval[1]*sr)
            tabDuration = interval[1] - interval[0]
            if (tabDuration < (valWin + valStep)):
                print("1: Segmento ID {:} de {:} seg é muito pequeno para processamento.".format(tabId,tabDuration))
                print("1: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
                
            tags = interval[2].split("-")
            # --- Trata digito errado
            if (tags[0].isdigit()):
                if (tags[0] == '1'):
                    tags[0] = 'i'
                elif (tags[0] == '3'):
                    tags[0] = 'e'
                elif (tags[0] == '4'):
                    tags[0] = 'a'
                else:
                    tags[0] = 'o'
            # ----------------------
            tags[0] = tags[0].replace(" ","")
            
            if not (len(tags) == 5):
                print("2: Etiqueta {:} ID {:} com problema no número de marcações.".format(interval[2],tabId))
                print("2: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                continue
            if (len(tags[2]) > 0):
                try:
                    nSib = int(tags[2])
                    if (nSib > maxNSyllab):
                        maxNSyllab = nSib
                except:
                    print("9: Problema com a TAG 2 como {:} de ID {:} ".format(tags[2],tabId))
                    continue
            else:
                continue
                    
            tabFonetica = tags[0]
            if (len(tabFonetica) == 1) and (tabFonetica.lower() in listConsoante):
                print("2: Etiqueta {:} ID {:} apresenta tag que nao e vogal {:}.".format(interval[2],tabId,tabFonetica))
                continue
            
            tabPalavra = tags[1]
            tabTonicidade = tags[3]
            tabDitongo = is_ditongo(tags[0])
                            
            selAudio = audio[nIni:nFim]
            form2, _ = format_lpc(selAudio,sr,winstep=valStep,winlen=valWin)
            inten = intensity(selAudio,sr,winstep=valStep,winlen=valWin)
            try:
                tabMeanHNR, _ = get_HNR(selAudio,sr,time_step=valStep,periods_per_window = 1.875)
            except:
                tabMeanHNR = 0;
                continue
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
            try:
                sibPosition = int(tags[4])
                if (int(tags[4]) >  len(grafPalavra)):
                    print("5: Etiqueta {:} ID {:} com problema. Possicao maior que o numero de silabas".format(interval[2],tabId))
                    print("5: Arquivo {:} no intervalo entre {:2.3f} e {:2.3f} segundos.".format(Path(tgFile).name,interval[0],interval[1]))
                    continue
            except:
                print("9: Problema com a TAG 4 como {:} de ID {:} ".format(tags[4],tabId))
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
            #phonTrancrib = ''.join(phonPalavra).replace("ˈ","")
            # TODO:
                # Nao identifica algumas vogais tipo em " ʧĩ"
            
            hasNasal = False
            idxNasal = []
            if ('͂' in phonSilaba) or ('̃' in phonSilaba):
                hasNasal = True    
                oriphonSilaba = phonSilaba
                phonSilaba = phonSilaba.replace('͂','').replace('̃','')
                print("Tratar as nasais")
                
            
            if (not has_vogal(tabFonologico)):
                print("7: ID {:} com problema, {:} sem vogal.".format(tabId, tabFonologico))
                continue
            
            tabPrecedente = 'NA'
            tabSeguinte = 'NA'
            if (tabDitongo):
                vogLen = len(tags[0])
                pos = grapSilaba.find(tags[0])
                if (pos == -1):
                    pos, valT = find_pos_of_tag(grapSilaba,tags[0])
                elif (pos == -1):
                    print("13: Problema na detecçao da tag ID {:}".format(tabId))
                    continue
                if (pos == 0) and (len(phonSilaba) == vogLen):
                    tabPrecedente = 'NA'
                    tabSeguinte = 'NA'
                if (pos == 0) and (len(phonSilaba) > vogLen):
                    tabPrecedente = 'NA'
                    tabSeguinte = phonSilaba[1]
                if (pos > 0) and (len(phonSilaba) == (pos+vogLen)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = 'NA'
                if (pos > 0) and (len(phonSilaba) > (pos+vogLen)):
                    tabPrecedente = phonSilaba[pos-1]
                    tabSeguinte = phonSilaba[pos+1]
            else:
                if (hasNasal):
                    print('Depurando...')
                pos = pos_indicated_vowel(phonSilaba, tags[0], tabId)
                if (pos == -1):
                    pos, valT = find_pos_of_tag(grapSilaba,tags[0])
                elif (pos == -1):
                    print("13: Problema na detecçao da tag ID {:}".format(tabId))
                    continue
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
            
            
            
            tabLetraPre = 'NA'
            tabLetraSeg = 'NA'
            # if (tabId == 10):
            #     print('Depurando...')
            if (tabDitongo):
                vogLen = len(tags[0])
                pos = tabPalavra.find(tags[0])
                if (pos == -1):
                    pos, valT = find_pos_of_tag(tabPalavra,tags[0])
                elif (pos == -1):
                    print("13: Problema na detecçao da tag ID {:}".format(tabId))
                    continue
                if (pos == 0) and (len(tabPalavra) == vogLen):
                    tabLetraPre = 'NA'
                    tabLetraSeg = 'NA'
                if (pos == 0) and (len(tabPalavra) > vogLen):
                    tabLetraPre = 'NA'
                    tabLetraSeg = tabPalavra[pos+vogLen]
                if (pos > 0) and (len(tabPalavra) == (pos+vogLen)):
                    tabLetraPre = tabPalavra[pos-1]
                    tabLetraSeg = 'NA'
                if (pos > 0) and (len(tabPalavra) > (pos+vogLen)):
                    tabLetraPre = tabPalavra[pos-1]
                    tabLetraSeg = tabPalavra[pos+vogLen]
            else:
                pos = pos_vowel_in_word(grafPalavra, vogalPos, tags[0])
                if (pos == 0) and (len(grafPalavra) == 1):
                    tabLetraPre = 'NA'
                    tabLetraSeg = 'NA'
                if (pos == 0) and (len(tabPalavra) > 1):
                    tabLetraPre = 'NA'
                    tabLetraSeg = tabPalavra[pos+1]
                if (pos > 0) and (len(tabPalavra) == (pos+1)):
                    tabLetraPre = tabPalavra[pos-1]
                    tabLetraSeg = 'NA'
                if (pos > 0) and (len(tabPalavra) > (pos+1)):
                    tabLetraPre = tabPalavra[pos-1]
                    tabLetraSeg = tabPalavra[pos+1]
                    
            if (hasNasal):
                print("recolocar a nasal")
                phonSilaba = oriphonSilaba
                
            # TODO: Retirar redundancia de tabOral e hasNasal
            tabOral = int((not hasNasal))
            tabFechada = int((grapSilaba[-1] in listConsoante))
            tabData = (tabId,tabDuration, tabF1,tabF2,tabF1_b,tabF2_b,tabIntensity,
                       interval[2],vogalPos,
                       tabMeanHNR, tabFonetica,tabFonologico, tabDitongo, 
                       tabPalavra, tabTonicidade, tabPrecedente, tabSeguinte,tabFechada,
                       nSib, tabOral, tabLetraPre.lower(), tabLetraSeg.lower(),
                       tabSexo,tabFileName)
            
            strData = "{:}\n".format(tabData).replace("(","").replace(")","").replace(" ","").replace(",","\t")
            # sys.exit("Saida de depuraçao")
            csvLines.append(strData)
            tabId = tabId + 1
        # sys.exit("Saida de depuraçao - RODOU APENAS CAMADA 1!")
        
with open(CSVFILE, 'w') as file:
        file.writelines(csvLines)

# Processamento de dados de DE vogais para ICPhS

### Cálculo de características acústicas em fricativas

Cálculo de características acústicas em vogais

Minha humilde contribuição para o trabalho de Maria Cantoni e Thais Cristofaro

Data: 28/12/2022 (Lula já eleito)

Autor: Adelino P. Silva

mail: adelinocpp@yahoo.com, adelinocpp@gmail.com, adelino@ufmg.br

tel: +55 31 98801-3605

---
Rotina para a leitura de arquivos do tipo .TextgGrid (praat) e calculo de 
taxas espectrais como:

1 - HNR (de Signal_Analysis.features.signal)

2 - centro de gravidade espectral do tipo 1, 2 e 2/3, (vide link abaixo) e

3 - Formantes

4 - Etc...
~~3 - Low-frequency-to-total intensities ratio (LTF), vide GRADOVILLE (2021).~~

Centro de gravidade praat:
https://www.fon.hum.uva.nl/praat/manual/Spectrum__Get_centre_of_gravity___.html

Utiliza silabificação do Pɛtɾʊs (PhonEtic TRanscriber for User Support), com alterações, originalmente desenvolvido por:

[alessandrobokan](https://github.com/alessandrobokan) (Obrigado. Desculpa se te divulguei sem permissão.)

e-mail: alessandro.bokan@gmail.com

github Pɛtɾʊs PT-BR: https://github.com/alessandrobokan/PETRUS

---
## Como rodar (testado em Ubuntu 20.20 e Mint 21.3)

Construa o sistema de diretórios
```
principal (qualquer diretorio)
    |-Audios
        |-Arquivo_audio_01.wav
        |-Arquivo_audio_01.TextGrid
        |-Arquivo_audio_02.wav
        |-Arquivo_audio_02.TextGrid
        ...
    |-Processamento (conteudo deste repositório)
        |-g2p
        |-stress
        |-syllables
        |-utils
        |-P00_Compute_Vogal_Features_v0.py
        ...
```
Instale o pacote (ainda faço um arquivo re requeriments.txt)

```
$ pip install Signal_Analysis
```

Execute no terminal (no dieretório onde está o arquivo P00_Compute_Vogal_Features_v0.py)

```
$ python3 P00_Compute_Vogal_Features_v0.py
```

Aguarde gerar o arquivo "csvDataAudios.csv". No caminho aparecem muitas mensagens de falhas, mas a maioria dá certo.

### Observações:

A variável "AUDIO_FOLDER" indica o diretŕorio dos arquivos de áudio e "CSVFILE" o arquivo CSV de saída.

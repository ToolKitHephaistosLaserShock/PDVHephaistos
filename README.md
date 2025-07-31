# Free and easy tool for PDV simple analysis 

Thanks to Gabriel Prudhomme. 

Étude du nuage de particules éjectées sous choc : apports de la Vélocimétrie Hétérodyne. Mécanique des matériaux [physics.class-ph]. 

Ecole nationale supérieure d'arts et métiers

ENSAM, 2014.

Français. ⟨NNT : 2014ENAM0044⟩. ⟨tel-01165754⟩

https://pastel.hal.science/tel-01165754

Base signal for analysis :Tension(Time) .csv file

**class > PDV(Time,Tension,ChainResponse,PDVShift,PDVFactor,FName,ShotNumber)**

    - Raw data : Time Tension
    - Chains response (GHz)
    - PDVShift 
    - PDVFactor (m/s/Hz)
    - FName - File name of raw datas   
    - ShotNumber - Files directory
*************************************************************
**Directory Structure :**

**Working directory**

    - PDVWorking.py
    - PDVExtractSignalAndAnalysisBetaTest.py
    -- ShotNumber Directory
        - ShotNumber.csv #data sets file .csv Time(s),Tension(V)
        - ShotNumber.csvVelocity.csv # Velocity extracted file .cvs Time(ns),Velocity(m/s)
        - ShotNumber.csvVelocity.png # Velocity plot image .png
        - ShotNumber.csvRawDatas.png # RawDatas plot image .png
        - ShotNumber.csvSpectrogram.png # Last calculation for Spectrogram plot image .png
******************************

def DataLoad(self,LinesSuppressed) load data from .csv file Tension(Time). 
       
def PDVSetFrAcquisition(self) - Extract sample rate in GS/s

def SetPDVFFT(self) - calculate FFT of raw datas :  Tension and related Time

def SetSTFTPDV(self,nperseg) - Calculate STFT from raw data on number of point - nperseg

def SetVelocity() - calculate velocity m/ss

def PDVReport() - pdf report with all datas and graph for basic analysis, datas, FFFT, Spectrogram,baseline.  

**Working on Conda env with  :**

    - matplotlib
    - numpy
    - pyqt
    - pyqtwebengine
    - python 3.11
    - scipy
    - spyder
    - tk
    - pandas
    - reportlab
    
**ShotTest.zip contains an example to be tested**


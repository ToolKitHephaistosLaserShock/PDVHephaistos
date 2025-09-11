# Free and easy tool for PDV simple analysis 

![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/main/DataLoad.png).

Thanks to Gabriel Prudhomme. 

Étude du nuage de particules éjectées sous choc : apports de la Vélocimétrie Hétérodyne. Mécanique des matériaux [physics.class-ph]. 

Ecole nationale supérieure d'arts et métiers

ENSAM, 2014.

Français. ⟨NNT : 2014ENAM0044⟩. ⟨tel-01165754⟩

https://pastel.hal.science/tel-01165754

Base signal for analysis by STFT and Wavelet:Tension(Time) .csv file

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
**Main functions**

    - def DataLoad(self,LinesSuppressed) load data from .csv file Tension(Time). 
      
    - def PDVSetFrAcquisition(self) - Extract sample rate in GS/s

    - def SetPDVFFT(self) - Calculate FFT of raw datas :  Tension and related Time

    - def SetSTFTPDV(self,nperseg) - Calculate STFT from raw data on number of point - nperseg

    - def SetWavelet(self,WidthWavelet) - Calculate Wavelet from raw data on number of point - nperseg and function (Morelet MexHat)

    - def SetVelocity() - Calculate velocity m/ss

    - def PDVReport() - Generate pdf report with all datas and graph for basic analysis, datas, FFFT, Spectrogram,baseline.  
	
	- def ExtractVelocityNotebookAuto() - Extrate velocity profile automatically, based on spectrum maximum
	
	- def BaseLineDelete() - Delete baseline

**Automatic Velocity Profile Extraction :**

	For each time, extract maximum frequency in the spectrogram at this time
	
**Automatic BaseLine Management :**

![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/main/Capture_BaseLine_Delete.png).

	- Delete : extract spectrum at time previous any shock arrival, substract this spectrum at all time, ponderate by baseline ratio at each time
	
	- Reset : Two Spectrogram variables exist. The first one is the spectrogram calculated while the second is the displayed spectrogram. First remains untouched. Calculations (such as deleting the baseline) are only made on the second. When "adding" the baseline, reload the calculated spectrogram as the displayed one. No calculations.


**Working on Conda env with :**

    - matplotlib
    - numpy
    - pyqt
	- PyWavelets
    - pyqtwebengine
    - python 3.11
    - scipy
    - spyder
    - tk
    - pandas
    - reportlab
    
**ShotTest.zip contains an example to be tested**


![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/main/SFTinteractive.png "SFT Interactive").

# Free and easy tool for PDV simple analysis 

## Extraction window in Hz and s for low velocity/pendulum 

Citation : Berthe, L., & Delalande, R. (2025). ToolKitHephaistosLaserShock/PDVHephaistos: PDVHephaistos (Latest). Zenodo. https://doi.org/10.5281/zenodo.16684803 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16684803.svg)](https://doi.org/10.5281/zenodo.16684803)


![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/OrientedLowVelocityPendulum/DataLoad.png).

Thanks to Gabriel Prudhomme. Étude du nuage de particules éjectées sous choc : apports de la Vélocimétrie Hétérodyne. Mécanique des matériaux [physics.class-ph]. 
Ecole nationale supérieure d'arts et métiers,ENSAM, 2014.Français. ⟨NNT : 2014ENAM0044⟩. ⟨tel-01165754⟩ https://pastel.hal.science/tel-01165754

**Base signal file for analysis : Tension,Time  .csv file** 

Comma can be changed easily in the DataLoad function

You can suppress a number of headlines from the base signal file

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
		
    - def SetWaveletTransformPDV(self, WidthWavelet)- Continuous Wavelet Transform (CWT) using PyWavelets - WidthWavelet : number of scale (scales).
    
    - def SetVelocity() - calculate velocity m/ss

    - def SetVelocity() - Calculate velocity m/ss

    - def PDVReport() - Generate pdf report with all datas and graph for basic analysis, datas, FFFT, Spectrogram,baseline.  
	
	- def ExtractVelocityNotebookAuto() - Extrate velocity profile automatically, based on spectrum maximum
	
	- def BaseLineDelete() - Delete baseline

**Automatic Velocity Profile Extraction :**
![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/056bb88ebd103ba899c32df5100b4d9f7d6c3603/Figure_Exemple_ExtractionF.png).

	For each time step, maximum frequency is identified. Frequency at half-amplitude is also detected on both left and right of the maximum.
	
**Automatic BaseLine Management :**
![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/c90a0e8b6aa25a6592d646e8c44c0f1792cc5357/Capture_BaseLine_Delete.png).

	- Delete : extract spectrum at time previous any shock arrival, substract this spectrum at all time, ponderate by baseline ratio at each time
	
	- Reset : Two Spectrogram variables exist. The first one is the spectrogram calculated while the second is the displayed spectrogram. First remains untouched. Calculations (such as deleting the baseline) are only made on the second. When "adding" the baseline, reload the calculated spectrogram as the displayed one. No calculations.

** Help in-app:**
![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/2fb299543817dc3a313a5cb576fec5741a7e6ae4/Capture_Help.png).

Quick help is available within the app. Click on the "?" question marks as in the picture.

**Working on Conda env with :**

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

![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/OrientedLowVelocityPendulum/SFTinteractive.png)

![](https://github.com/ToolKitHephaistosLaserShock/PDVHephaistos/blob/OrientedLowVelocityPendulum/Velocity.png)

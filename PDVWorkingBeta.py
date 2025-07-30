#!/usr/bin/env python
# coding: utf-8
"""
Thanks to Gabriel Prudhomme. 
Étude du nuage de particules éjectées sous choc : apports de la Vélocimétrie Hétérodyne. Mécanique des matériaux [physics.class-ph]. 
Ecole nationale supérieure d'arts et métiers
ENSAM, 2014. 
Français. ⟨NNT : 2014ENAM0044⟩. ⟨tel-01165754⟩
https://pastel.hal.science/tel-01165754


Base signal for analysis :Tension(Time)
class > PDV(Time,Tension,ChainResponse,PDVShift,PDVFactor,FName,ShotNumber)
    -Raw data : Time Tension
    -Chains response (GHz)
    -PDVShift 
    -PDVFactor (m/s/Hz)
    -FName - File name of raw datas   
    -ShotNumber - Files directory
*************************************************************
Directory Structure : 
Working directory
    -PDVWorking.ipynb
    --ShotNumberDirectory
        -RawDatas
        -Graphs 
        -Report
******************************
def DataLoad(self,LinesSuppressed) load data from .csv file Tension(Time). 
       
def PDVSetFrAcquisition(self) - Extract sample rate in GS/s
def SetPDVFFT(self) - calculate FFT of raw datas :  Tension and related Time
def SetSTFTPDV(self,nperseg) - Calculate STFT from raw data on number of point - nperseg
def SetVelocity(s%matplotlib inlineelf) - calculate velocity m/ss
def PDVReport(self) - pdf report with all datas and graph for basic analysis, datas, FFFT, Spectrogram, baseline.  
"""
from PDVExtractSignalAndAnalysisBetaTest import *

C=3e8 #m/s


## Input PDV parameters 
LambdaLaser=1550 #nm
LambdaReference=1550 #nm
#PDV Shift
PDVShift=1#GHz

# Reponse chain
ChainResponse=1#Ghz

nperseg =2000 # FFT glissante

#File Name
FName = "ShotTest"
ShotNumber='ShotTest' #Directory Shot

#SI
LambdaLaser=LambdaLaser*1e-9
FrLaser= C/LambdaLaser
LambdaReference=LambdaReference*1e-9
FrReference= C/LambdaReference
print("Lambda Laser : (m)     :",f"{LambdaLaser}")
print("FrLaser : (Hz)         :",f"{FrLaser:e}")
print("Lambda Reference : (m) :",f"{LambdaReference}")
print("FrReference : (Hz)     :",f"{FrReference:e}")
PDVFactor=LambdaLaser/2
print("PDVFactor (m)          :",f"{PDVFactor:e}")



#Extraction from .trc generate .csv datas set (Scope signal from lecroy)
#trc=Trc()
#trc.open(FName)
#print("Base Time : ", trc.ScopeStatus["TIMEBASE"])
#print ("Number of points :",trc.ScopeStatus["WAVE_ARRAY_COUNT"])

#Open class for ShotNumber
ShotNumber=PDV(ChainResponse,PDVShift,PDVFactor,FName,ShotNumber,nperseg)
ShotNumber.NotebookGraph()
ShotNumber.runSTFTPDVInteractive()

###Directory Shots
#print ("##Goto directory ShotNUmber")
#WorkDirectory=ShotNumber
#print('Ask WorkDirectory : ',WorkDirectory)
#print('Current Directory before change : ', os.getcwd())
#os.chdir(WorkDirectory)
#print('Current Directory : ', os.getcwd())

#### Basic analysis > in report

##ShotNumber.DataLoad(1)
#calculate Acquision time scale
#ShotNumber.PDVSetFrAcquisition()
#FFT                 
#ShotNumber.SetPDVFFT()
# Calculate STFT
#ShotNumber.SetSTFTPDV(nperseg)
#Calculate velocity
#ShotNumber.SetVelocity()
#graph spectrogram & velocty
#ShotNumber.GraphSpectrogram()
#ShotNumber.WindowDatas()
#ShotNumber.CreateSTFTPDVInteractive()
#ShotNumber.runSTFTPDVInteractive()
#ShotNumber.CreateSTFTPDVInteractive ()
#ShotNumber.runSTFTPDVInteractive()
#Extraction manualy point point point
#ShotNumber.ExtractVelocityProfile(ShotNumber.Time_stft,ShotNumber.Velocity,ShotNumber.PDVSpectrogram)
#ShotNumber.ExtractVelocityProfileGraphSave()
#ShotNumber.FrequenceDominant(1,ShotNumber.PolygonTime,ShotNumber.PolygonVelocity,ShotNumber.PolygonPDVSpectrogram)
#ShotNumber.FillBetweenVelocityProfile(20,20)
#ShotNumber.PDVReport()








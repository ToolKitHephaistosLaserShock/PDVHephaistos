#!/usr/bin/env python
# coding: utf-8

import numpy as np

import pywt

from scipy.signal import stft

from pylab import *  # ? used for datetime ?

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

import csv

import time

import os
import sys

import tkinter.font as tkFont
from tkinter import ttk
import tkinter as tk
import tkinter.filedialog as fd

import struct

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer,KeepTogether,tables,PageBreak
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.pagesizes import A4,landscape
from reportlab.lib.units import inch,cm,mm
from reportlab.pdfgen import canvas

C=3e8 #m/s

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
    -PDVWorking.py
    --ShotNumberDirectory
        -RawDatas
        -Graphs 
        -Report
******************************
def DataLoad(self,LinesSuppressed) load data from .csv file Tension(Time). 
def PDVSetFrAcquisition(self) - Extract sample rate in GS/s
def SetPDVFFT(self) - calculate FFT of raw datas :  Tension and related Time

def SetSTFTPDV(self,nperseg) - Calculate STFT from raw data on number of point - nperseg
def SetVelocity(s%matplotlib inlineelf) - calculate velocity m/s
def PDVReport(self) - pdf report with all datas and graph for basic analysis, datas, FFFT, Spectrogram, baseline.  
"""

class RedirectConsole:
    
    def __init__(self, text_widget):
        self.output = text_widget
        
    def write(self, string):
        self.output.insert(tk.END, string)
        self.output.see(tk.END)
        
    def flush(self):
        pass
    

class PDV :
    def __init__(self,LambdaLaser,ChainResponse,Shift,FName,ShotNumber,nperseg,WidthWavelet):
        
        #all are calculate in SI but on application print in Nm and Ghz
        self.LambdaLaser=LambdaLaser*1e-9 #>SI
        self.Shift=Shift*1e9 #Hz PDV shift of reference Laser
        self.FName=FName #ShotName
        self.ChainResponse=ChainResponse*1e9
        self.ShotNumber=ShotNumber
        self.PDVFactor=self.LambdaLaser/2
        self.nperseg=nperseg
        self.WidthWavelet=WidthWavelet
        self.STFTPDVWindow='hamming'
        self.WaveletFunctionPDV='morl'
				
        self.nperseg=500
        
			 	# Tools functions for calculation *****************************
        #print ("Design*****")
        self.VPivot=self.LambdaLaser*self.Shift/2
        #print("VPivot (m/s) :",self.VPivot)
        self.MaxVelocityForChainResponse=self.PDVFactor*self.ChainResponse
        #Close all fig
        plt.close('all')
        
    
    def DataLoad(self,LinesSuppressed):
        print("## Data Loading From .csv with , as separator ###############")
        print("Suppressed lines:", LinesSuppressed)
    
        self.Time = []  # Use lists initially
        self.Tension = [] 
        name=self.FName
        print (name)
        DataSet=csv.reader(open(name),delimiter=',')
        
         # Ouvrir le fichier avec 'with' pour garantir la fermeture automatique
        for i,e in enumerate(DataSet): 
            if i>LinesSuppressed : 
                ti,vi=float(e[0]),float(e[1])
                self.Time=np.append(self.Time,ti) ####ICICICICICICI
                self.Tension=np.append(self.Tension,vi)
                
        print("Number of points from DataLoad:", len(self.Time))
        return
    
    def PDVSetFrAcquisition(self):
        #Data extraction of acquisition sample rate in Sample/s
        self.Dtime = (self.Time[-1] - self.Time[0])/(len(self.Time) - 1)
        
        self.Lbl_Time_unit_Long =  self.Cbo_TimeUnits.get()
        
        # self.Choices_time_unit = ['ns', 'µs', 'ms', 's']      # Possible unit to select
        
        if self.Lbl_Time_unit_Long == self.Choices_time_unit[0]: # ns
            self.Time_OM_Long = 1e-9
            if (self.Dtime/self.Time_OM_Long)<1E-3:
                self.Lbl_Time_unit_Court = 'ps'
                self.Time_OM_Court = 1e-12
                self.Lbl_Freq_unit = 'THz'
            else:
                self.Lbl_Time_unit_Court = self.Lbl_Time_unit_Long
                self.Time_OM_Court = self.Time_OM_Long
                self.Lbl_Freq_unit = 'GHz'
        elif self.Lbl_Time_unit_Long == self.Choices_time_unit[1]: # µs
            self.Time_OM_Long = 1e-6
            if (self.Dtime/self.Time_OM_Long)<1E-3:
                self.Lbl_Time_unit_Court = 'ns'
                self.Time_OM_Court = 1e-9
                self.Lbl_Freq_unit = 'GHz'
            else:
                self.Lbl_Time_unit_Court = self.Lbl_Time_unit_Long
                self.Time_OM_Court = self.Time_OM_Long
                self.Lbl_Freq_unit = 'MHz'
        elif self.Lbl_Time_unit_Long == self.Choices_time_unit[2]: # ms
            self.Time_OM_Long = 1e-3
            if (self.Dtime/self.Time_OM_Long)<1E-3:
                self.Lbl_Time_unit_Court = 'µs'
                self.Time_OM_Court = 1e-6
                self.Lbl_Freq_unit = 'MHz'
            else:
                self.Lbl_Time_unit_Court = self.Lbl_Time_unit_Long
                self.Time_OM_Court = self.Time_OM_Long
                self.Lbl_Freq_unit = 'kHz'
        elif self.Lbl_Time_unit_Long == self.Choices_time_unit[3]: # s
            self.Time_OM_Long = 1
            if (self.Dtime/self.Time_OM_Long)<1E-3:
                self.Lbl_Time_unit_Court = 'ms'
                self.Time_OM_Court = 1e-3
                self.Lbl_Freq_unit = 'kHz'
            else:
                self.Lbl_Time_unit_Court = self.Lbl_Time_unit_Long
                self.Time_OM_Court = self.Time_OM_Long
                self.Lbl_Freq_unit = 'Hz'
                
        self.Freq_OM = round(1/self.Time_OM_Court)
        
        print("## PDVSetFrAcquisition calculation")
        
        self.FAcquisition = 1/ self.Dtime
        print("Dtime (assumed in s) :" + str(self.Dtime))
        print("FAcquisition (assumed in Hz) :" + str( self.FAcquisition))
        
        print("Lbl_Time_unit_Long : " + str(self.Lbl_Time_unit_Long))
        print("self.Time_OM_Long : " + str(self.Time_OM_Long))
        print("(self.Dtime/self.Time_OM_Long) : " + str((self.Dtime/self.Time_OM_Long)))
        print("self.Lbl_Time_unit_Court : " + str(self.Lbl_Time_unit_Court))
        print("self.Time_OM_Court : " + str(self.Time_OM_Court))
        print("self.Lbl_Freq_unit : " + str(self.Lbl_Freq_unit))
        print("self.Freq_OM :" + str(self.Freq_OM))

        print ("## PDVSetFrAcquisition calculation")
        self.Dtime=self.Time[2]-self.Time[1]
        self.FAcquisition=1/self.Dtime
        print ("DTime (ns) :", self.Dtime*1e9)
        print ("FAquisition (GS/s) :",f"{self.FAcquisition*1e-9:e}")
        return
        
    
    def SetPDVFFT(self):
        #Calculate simple FFT raw data
        print ("## FTT signal calculation")
        self.N = len(self.Time)  # Taille du signal
        print ("Number of points: ", self.N)
        self.HSignalFFT= np.fft.rfft(self.Tension)  # Calcul de la FFT
        self.PDVSignalFFTTime= np.fft.rfftfreq(self.N, 1/self.FAcquisition)  # Axe fréquentiel
        return
    
    def SetSTFTPDV(self,nperseg):
        self.WindowsSize=self.Dtime*self.nperseg
        self.FePDV, self.Time_stft, self.PDVSpectrogram = stft(self.Tension, self.FAcquisition, nperseg=self.nperseg,window=self.STFTPDVWindow)
        return
    
    def SetWaveletTransformPDV(self, WidthWavelet):
        """
        Continuous Wavelet Transform (CWT) using PyWavelets.
        WidthWavelet : number of scale (scales).
        """
        print ("## Wavelet signal calculation")
        # Define scale
        scales = np.arange(1, WidthWavelet)
        # Calcul CWT 
        self.WaveletFrequencies=[]
        self.WaveletSignalPDV=[]
        self.WaveletSignalPDV, self.WaveletFrequencies = pywt.cwt(
            self.Tension, scales, self.WaveletFunctionPDV, sampling_period=self.Dtime
        )
        return
        
    
    def SetVelocity(self):
        #Calculate PDVFactor
        print ("## Velocity signal calculation")
        self.Velocity=self.PDVFactor*self.FePDV
        return
    
 # Tools functions for boxes and operations *****************************    
    def NotebookGraph(self):   
        self.root = tk.Tk()
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(family="Verdana", size=12)
        self.root.title("PDV Analysis")
        self.shot_var = tk.StringVar()
        self.fname_var= tk.StringVar()
        self.nperseg_var=tk.StringVar()
        self.WidthWavelet_var=tk.StringVar()
        self.shot_dir=tk.StringVar()    
        self.ChainResponse_var=tk.StringVar()
        
        # Screen size
        self.Screen_H = self.root.winfo_screenheight()
        self.Screen_W = self.root.winfo_screenwidth()
        self.px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        print("\n width x height = %d x %d (in pixels)\n" %(self.Screen_W, self.Screen_H))
        
        # Onglets
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=1)
        
        # Onglet Inputs/inputs - design PDV 
        self.frame_inputs = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_inputs, text="Datas Load & Operations")
        self.CreateInputsTab(self.frame_inputs)
        
        #Console output
        self.CreateConsoleTab()
        sys.stdout = RedirectConsole(self.text_console)
        sys.stderr = RedirectConsole(self.text_console)
        
        # Close all figures
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
						
        
    def CreateInputsTab(self, parent): # Tab for data set input and PDV calculation
        style = ttk.Style()
        style.configure('TButton', font=('Arial', 12, 'bold'))
        
        # Variables for tab printing on screen
        self.LineSuppressed_var= tk.StringVar(value=1)
        self.shot_dir = tk.StringVar()
        self.fname = tk.StringVar()
        self.nperseg_var = tk.IntVar(value=self.nperseg)
        self.ChainResponse_var = tk.DoubleVar(value=self.ChainResponse * 1e-9)
        #Unit for print in nm
        Wavelength=self.LambdaLaser*1e9
        self.LambdaLaser_var= tk.StringVar(value=f'{Wavelength:.3f}')
        self.Shift_var= tk.DoubleVar(value=self.Shift * 1e-9)
        self.MaxVelocityForChainResponse_var=tk.StringVar(value=f'{self.MaxVelocityForChainResponse:.3f}')
        self.VPivot_var=tk.StringVar(value=f'{self.VPivot:.3f}')
        
        ttk.Label(parent, text="PDV Analysis", font=("Arial", 14, "bold")).pack(anchor="w", padx=10)
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=2, ipady=3)
        #Datas selection ********************************************************
        # Directory choice (ShotNumber)
        ttk.Label(parent, text="Shot Directory:").pack(pady=5)
        frame_dir = tk.Frame(parent)
        frame_dir.pack()
        tk.Entry(frame_dir, textvariable=self.shot_dir, width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame_dir, text="Select Directory", command=self.select_directory).pack(side=tk.LEFT)
    
        # File date choice (FName)
        ttk.Label(parent, text="File to analyse .csv with template (Time(s),Tension(V))").pack(pady=2)
        frame_file = tk.Frame(parent)
        frame_file.pack()
        tk.Entry(frame_file, textvariable=self.fname, width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame_file, text="Select File ", command=self.select_file).pack(side=tk.LEFT)
        
        # Unit (time) selection
        self.Lbl_Time_Unitstmp = 'ns'
        self.List_Lbl_Time_Units = tk.StringVar(value=self.Lbl_Time_Unitstmp)
        self.Choices_time_unit = ['ns', 'µs', 'ms', 's']      # Possible unit to select
        
        self.Frame_Time_Unit = tk.Frame(parent)
        self.Frame_Time_Unit.pack(side=tk.TOP)
        
        self.BtnHelp_Time_Unit = tk.Button(self.Frame_Time_Unit, text="?", command=self.InterfaceHelpTimeUnit)
        self.BtnHelp_Time_Unit.pack(anchor="w", side=tk.LEFT)
        ttk.Label(self.Frame_Time_Unit, text="Time unit to work with:").pack(side=tk.LEFT, anchor="w")
        self.Cbo_TimeUnits = ttk.Combobox(self.Frame_Time_Unit, textvariable=self.List_Lbl_Time_Units, values=self.Choices_time_unit, state="readonly")
        self.Cbo_TimeUnits.pack(anchor = "w", padx=5, side=tk.LEFT)
        ttk.Label(self.Frame_Time_Unit, text="assuming file data are in s.").pack(side=tk.LEFT, anchor="w")

        ttk.Label(parent, text="Lines Suppressed on data file").pack(pady=5)
        frame_LineSuppressed= tk.Frame(parent)
        frame_LineSuppressed.pack()
        tk.Entry(parent, textvariable= self.LineSuppressed_var, width=15).pack()
        
        # Launch analysis
        ttk.Label(parent, text="Spectrogram, Raw datas and velocity figures are saved in png format ").pack(pady=10)
        ttk.Label(parent, text="Velocity data set in .csv file in ShotNumber Directory").pack(pady=5)
        ttk.Button(parent, text="Load Data Set for analysis", style='TButton',command=self.launch_analysis).pack(pady=5)
        
        ttk.Label(parent, text="PDV parameters", font=("Arial", 14, "bold")).pack(anchor="w", padx=10)
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=2)
        
        # LaserPDV Wavelength
        ttk.Label(parent, text="Laser PDV Wavelength (nm) :").pack(pady=10)
        tk.Entry(parent, textvariable=self.LambdaLaser_var, width=15).pack()
        
        # Line for tab
        line_frame = ttk.Frame(parent)
        line_frame.pack(anchor="w", pady=5)
        
        # Bloc 1 : Chain Response
        chain_frame = ttk.Frame(line_frame)
        chain_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(chain_frame, text="Chain Response (Ghz) :").pack(anchor="w")
        tk.Entry(chain_frame, textvariable=self.ChainResponse_var, width=15).pack()
        
        # Bloc 1 : Max Velocity
        velocity_frame = ttk.Frame(line_frame)
        velocity_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(velocity_frame, text="Max Velocity (m/s)").pack(anchor="w")
        tk.Entry(velocity_frame, text=self.MaxVelocityForChainResponse_var, width=15).pack()
        
        line_frame = ttk.Frame(parent)
        line_frame.pack(anchor="w", pady=5)
        
        # Bloc 1 : Shift
        chain_frame = ttk.Frame(line_frame)
        chain_frame.pack(side=tk.LEFT, padx=10)
        
        # Laser Shift
        ttk.Label(chain_frame, text="LaserPDV Shift in Ghz :").pack(anchor="w")
        tk.Entry(chain_frame, textvariable=self.Shift_var, width=15).pack()
        
        # Bloc 2 : Pivot Velocity
        velocityPivot_frame = ttk.Frame(line_frame)
        velocityPivot_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(velocityPivot_frame, text="Pivot Velocity (m/s)").pack(anchor="w")
        tk.Entry(velocityPivot_frame, text=self.VPivot_var, width=15).pack()
        
        # PDV parameters lauch
        
        ttk.Button(parent, text="PDV Parameters Calculations", style='TButton',command=self.PDVParameters).pack(pady=10)
        
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=5, ipady=3)
        ttk.Button(parent, text="Exit", style='TButton', command=lambda: os._exit(0)).pack(pady=5)
        
    
    def InterfaceHelpTimeUnit(self):
        self.Inter_Help_TimeUnit = tk.Tk()
        
        self.LblTitle_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="Define time unit", font=("Arial", 14, "bold"))
        self.LblTitle_Help_TimeUnit.pack(anchor="w", pady=5)
        
        self.Lbl1_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="Select the order of magnitude corresponding to the duration of the signal.")
        self.Lbl1_Help_TimeUnit.pack(anchor="w")
        
        self.Lbl2_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="Interface information will be given in this unit.")
        self.Lbl2_Help_TimeUnit.pack(anchor="w")
        
        self.Lbl3_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="Data are assumed to be given in s in the selected file.")
        self.Lbl3_Help_TimeUnit.pack(anchor="w")
        
        self.Lbl4_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="If time step is detected to be below 1/1000th of the selected unit, some indications will be given with lower order of magnitude.")
        self.Lbl4_Help_TimeUnit.pack(anchor="w")
        
        self.Lbl5_Help_TimeUnit = tk.Label(self.Inter_Help_TimeUnit, text="Frequencies such as acquisition are given accordingly.")
        self.Lbl5_Help_TimeUnit.pack(anchor="w")
        
        self.Btn_QuitHelp_ManVelExtract = tk.Button(self.Inter_Help_TimeUnit, text="Leave", command=lambda: self.Inter_Help_TimeUnit.destroy())
        self.Btn_QuitHelp_ManVelExtract.pack(anchor="w")
        
        self.Inter_Help_ManVelExtract.mainloop()
    def select_directory(self):
        dirname = fd.askdirectory(title="Select Shot Directory")
        if dirname:
            self.shot_dir.set(dirname)
        
    
    def select_file(self):
        initialdir = self.shot_dir.get() if self.shot_dir.get() else "."
        filename = fd.askopenfilename(title="Select raw datas file .csv (Time(s), Tension(V))", initialdir=initialdir)
        if filename:
            self.fname.set(os.path.basename(filename))
            self.selected_file_fullpath = filename  # Full path if mandatory later
				
		
    def CreateConsoleTab(self):
        self.frame_console = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_console, text="Console")
        self.text_console = tk.Text(self.frame_console, height=15, width=80)
        self.text_console.pack(fill='both', expand=True)
        
        sys.stdout = RedirectConsole(self.text_console)
        sys.stderr = RedirectConsole(self.text_console)
        
    
    def launch_analysis(self):
        
        for tab_id in self.notebook.tabs():
            tab_text = self.notebook.tab(tab_id, "text")
            if tab_text not in ("Datas Load & Operations","Console"):
                self.notebook.forget(tab_id)
        
        try:
        # Read width
            width_str = self.LineSuppressed_var.get()
            self.LinesSuppressed = int(float(width_str))
            if self.LinesSuppressed < 1:
                print("Lines Suppressed need to be > 0.")
                return
        except ValueError:
            print("Invalid value")
            return

        self.ShotNumber = self.shot_dir.get()
        self.FName = self.fname.get()
        self.nperseg = self.nperseg_var.get()
				
        #data print on consol output
        ##Directory Shots
        print ("##Goto directory ShotNUmber")
        WorkDirectory=self.ShotNumber
        print(self.ShotNumber,self.FName)
        print('Ask WorkDirectory : ',WorkDirectory)
        print('Current Directory before change : ', os.getcwd())
        os.chdir(WorkDirectory)
        print('Current Directory : ', os.getcwd())
        
        print ("## Start analysis#####")
        print ("Reminder parameter")
        print("Shot Number        : ",self.ShotNumber)
        print("FName :            : ",self.FName)
        print("ChainResponse (Hz) : ",self.ChainResponse,' Max Corresponding Velocity  (m/s) : ',self.ChainResponse*self.PDVFactor)
        print("Shift (Hz)         : ",self.Shift, ' Max Corresponding Velocity (m/s) : ',self.Shift*self.PDVFactor)
        print("PDVFactor m/s/Hz   : ",self.PDVFactor)
        print("window type        : " +self.STFTPDVWindow)
        print("Wavelet function   : " +self.WaveletFunctionPDV)
        print("Wavelet width      : ", self.WidthWavelet)
        print( "Lines suppressed  : ", self.LinesSuppressed)
        #Get datas and inital calculation
        self.DataLoad(self.LinesSuppressed)  # data load
				
        self.PDVSetFrAcquisition()
        self.SetPDVFFT()
        print ("self.nperseg      : ", self.nperseg, self.Dtime)
        self.SetSTFTPDV(self.nperseg)
        
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
    
        # tab Raw Datas
        self.frame_graphs = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_graphs, text="Raw Datas")
        self.NotebookGraphSpectrogram(self.frame_graphs)
    
        # tab STFT Interactive
        self.frame_stft = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_stft, text="STFT")
        self.CreateSTFTPDVInteractive(self.frame_stft)
        
        # tab Wavelet Interactive
        self.frame_Wavelet = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_Wavelet, text="Wavelet")
        self.CreateWaveletPDVInteractive(self.frame_Wavelet)
        
        # on tab STFT first
        self.notebook.select(self.frame_stft)
        
        #clean figures
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
						
    #interactive wavelet calculation analysis in playing with width and and function
    def CreateWaveletPDVInteractive(self, parent):
    
        self.wfig, self.wx = plt.subplots(figsize=(3, 2))
        self.wcanvas = FigureCanvasTkAgg(self.wfig, master=parent)
        self.wcanvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.wtoolbar = NavigationToolbar2Tk(self.wcanvas, parent)
        self.wtoolbar.update()
        self.wtoolbar.pack(side=tk.TOP, fill=tk.X)
    
        # self.wlabel1 = tk.Label(parent, text=f"Wavelet Window (Width) : {self.WidthWavelet} pt - Function : {self.WaveletFunctionPDV}")
        self.wlabel1 = tk.Label(parent, text=f"Wavelet Window (Width) : {self.WidthWavelet} pt - Function : {self.WaveletFunctionPDV}")
        self.wlabel1.pack()
        
        # self.wlabel2 = tk.Label(parent, text=f"Number of points : {len(self.Time)} pt, FAcquisition (GS/s): {self.FAcquisition*1e-9:e}")
        self.wlabel2 = tk.Label(parent, text=f"Number of points : {len(self.Time)} pt, FAcquisition ("+str(self.Lbl_Freq_unit)+"): "+str(self.FAcquisition/self.Freq_OM))
        # text=f"Number of points : {len(self.Time)} pt, FAcquisition ("+str(self.Lbl_Freq_unit)+"): "+str(self.FAcquisition/self.Freq_OM)
        self.wlabel2.pack()
        
        self.wx.set_title("Spectrogram " + self.FName)
        # self.wx.set_xlabel("Time (s)") #ICITIME
        self.wx.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        self.wx.set_ylabel("Frequency ("+str(self.Lbl_Freq_unit)+")")
        # self.wx.set_ylabel("Frequency (Hz)") #ICITIME
        
        self.LblTitle_ParamWL = tk.Label(parent, text="WL Parameters", font=("Arial", 14, "bold"))
        self.LblTitle_ParamWL.pack(anchor = "w", pady=5)
        
        self.Frame_ParamWL = ttk.Frame(parent)
        self.Frame_ParamWL.pack(anchor="w")
        
        self.WaveletFunctionPDV_var = tk.StringVar(value=self.WaveletFunctionPDV)
        Functions = ['morl', 'mexh']
        
        ttk.Label(self.Frame_ParamWL, text="Functions :").pack(anchor="w", side=tk.LEFT)
        combo = ttk.Combobox(self.Frame_ParamWL, textvariable=self.WaveletFunctionPDV_var, values=Functions, state="readonly")
        combo.pack(anchor="w", side=tk.LEFT)
        combo.bind('<<ComboboxSelected>>', lambda e: self.update_WaveletDVInteractiveplot(self.width_entry.get()))
        
        self.width_var = tk.StringVar(value=str(self.WidthWavelet))
        ttk.Label(self.Frame_ParamWL, text="Wavelet width :").pack(anchor="w", side=tk.LEFT)
        
        self.width_entry = ttk.Entry(self.Frame_ParamWL, textvariable=self.width_var, width=10)
        self.width_entry.pack(fill=tk.X, side=tk.LEFT)
        self.width_entry.bind("<Return>", self.update_WaveletDVInteractiveplot)
        
        ### Manual extraction wavelet
        self.Frame_ManualExtractWL1 = ttk.Frame(parent)
        self.Frame_ManualExtractWL1.pack(anchor="w", pady=5)
        
        self.BtnHelp_ManualExtractWL = tk.Button(self.Frame_ManualExtractWL1, text="?", command=self.InterfaceHelpManVelExtr)
        self.BtnHelp_ManualExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.LblTitle_ManualExtractWL = tk.Label(self.Frame_ManualExtractWL1, text="Manual Velocity Extraction", font=("Arial", 14, "bold"))
        self.LblTitle_ManualExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.Frame_ManualExtractWL2 = ttk.Frame(parent)
        self.Frame_ManualExtractWL2.pack(anchor="w")
        
        self.Instruction_ManualExtractWL = tk.Label(self.Frame_ManualExtractWL2, text="Manually extract velocity by clicking on the spectrogram.")
        self.Instruction_ManualExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.Btn_ManualExtractWL = tk.Button(self.Frame_ManualExtractWL2, text="Start Extraction", command=self.ExtractWaveletVelocityNotebook)
        self.Btn_ManualExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        ### Automatic extraction wavelet
        self.Frame_AutoExtractWL1 = ttk.Frame(parent)
        self.Frame_AutoExtractWL1.pack(anchor="w", pady=5)
        
        self.BtnHelp_AutoExtractWL = tk.Button(self.Frame_AutoExtractWL1, text="?", command=self.InterfaceHelpManVelExtr)
        self.BtnHelp_AutoExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.LblTitle_AutoExtractWL = tk.Label(self.Frame_AutoExtractWL1, text="Automatic Velocity Extraction", font=("Arial", 14, "bold"))
        self.LblTitle_AutoExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.Frame_AutoExtractWL2 = ttk.Frame(parent)
        self.Frame_AutoExtractWL2.pack(anchor="w")
        
        self.Instruction_AutoExtractWL = tk.Label(self.Frame_AutoExtractWL2, text="Automatic extraction in a defined window.")
        self.Instruction_AutoExtractWL.pack(anchor="w", side=tk.LEFT, pady=5)
        
        self.Frame_MinFreqAutoVelExtractWL = ttk.Frame(self.Frame_AutoExtractWL2)
        self.Frame_MinFreqAutoVelExtractWL.pack(anchor="w", side=tk.LEFT, padx=5)
        ttk.Label(self.Frame_MinFreqAutoVelExtractWL, text="Min Freq ("+str(self.Lbl_Freq_unit)+")").pack(anchor="w", side=tk.TOP)
        self.EntMinFreqWL = tk.Entry(self.Frame_MinFreqAutoVelExtractWL, width=15)
        self.EntMinFreqWL.pack(anchor="w", side=tk.TOP)
        
        self.Frame_MaxFreqAutoExtractWL = ttk.Frame(self.Frame_AutoExtractWL2)
        self.Frame_MaxFreqAutoExtractWL.pack(anchor="w", side=tk.LEFT, padx=5)
        ttk.Label(self.Frame_MaxFreqAutoExtractWL, text="Max Freq ("+str(self.Lbl_Freq_unit)+")").pack(anchor="w", side=tk.TOP)
        self.EntMaxFreqWL = tk.Entry(self.Frame_MaxFreqAutoExtractWL, width=15)
        self.EntMaxFreqWL.pack(side=tk.TOP)
        
        self.Frame_MinTimeAutoVelExtractWL = ttk.Frame(self.Frame_AutoExtractWL2)
        self.Frame_MinTimeAutoVelExtractWL.pack(anchor="w", side=tk.LEFT, padx=5)
        ttk.Label(self.Frame_MinTimeAutoVelExtractWL, text="T min ("+str(self.Lbl_Time_unit_Long)+")").pack(anchor="w", side=tk.TOP)
        self.EntMinTimeWL = tk.Entry(self.Frame_MinTimeAutoVelExtractWL, width=15)
        self.EntMinTimeWL.pack(side=tk.TOP)
        
        self.Frame_MaxTimeAutoVelExtractWL = ttk.Frame(self.Frame_AutoExtractWL2)
        self.Frame_MaxTimeAutoVelExtractWL.pack(anchor="w", side=tk.LEFT, padx=5)
        ttk.Label(self.Frame_MaxTimeAutoVelExtractWL, text="T max ("+str(self.Lbl_Time_unit_Long)+")").pack(anchor="w", side=tk.TOP)
        self.EntMaxTimeWL = tk.Entry(self.Frame_MaxTimeAutoVelExtractWL, width=15)
        self.EntMaxTimeWL.pack(side=tk.TOP)
        
        self.Btn_AutoExtractWL = tk.Button(self.Frame_AutoExtractWL2, text="Automatic extraction", command=self.AutoExtractVelocity_WL)
        self.Btn_AutoExtractWL.pack(anchor="w", side=tk.LEFT, padx=5)
        
        self.Lbl_AutoExtract_ErrorLblWL = ttk.Label(self.Frame_AutoExtractWL2, text="", foreground="red")
        self.Lbl_AutoExtract_ErrorLblWL.pack(anchor="w", side=tk.LEFT, padx=5)
        
        # Initial calculation
        self.SetWaveletTransformPDV(self.WidthWavelet)
        extent = [self.Time.min()/self.Time_OM_Long, self.Time.max()/self.Time_OM_Long,
              self.WaveletFrequencies.min()/self.Freq_OM, self.WaveletFrequencies.max()/self.Freq_OM]
        
        # Plot initial state
        plt.colorbar(self.wx.imshow(np.abs(self.WaveletSignalPDV),
               extent=extent,
               cmap='PRGn',
               aspect='auto',
               vmax=abs(self.WaveletSignalPDV).max()
               ))
        self.wcanvas.draw_idle()
			
    def update_WaveletDVInteractiveplot(self,event=None):
        try:
        # Read width
            width_str = self.width_var.get()
            self.WidthWavelet = int(float(width_str))
            if self.WidthWavelet <= 0:
                print("WidthWavelet need to be > 0.")
                return
        except ValueError:
            print("Invalid value")
            return        
        self.WaveletFunctionPDV = self.WaveletFunctionPDV_var.get()
        
        # Up date value in comments
        self.wlabel1.config(
        text=f"Wavelet Window (Width) = {self.WidthWavelet} pt - Function : {self.WaveletFunctionPDV}"
        )
        self.wlabel2.config(
        text=f"Number of points : {len(self.Time)} pt, FAcquisition ("+str(self.Lbl_Freq_unit)+"): "+str(self.FAcquisition/self.Freq_OM)
        )
        wxlim = self.wx.get_xlim()
        wylim = self.wx.get_ylim()

        # Update analysis with new value 
        self.SetWaveletTransformPDV(self.WidthWavelet)

        # Clear figures
        self.wx.clear()

        # Aplat time/frequency
        extent = [self.Time.min()/self.Time_OM_Long, self.Time.max()/self.Time_OM_Long,
              self.WaveletFrequencies.min()/self.Freq_OM, self.WaveletFrequencies.max()/self.Freq_OM]

        self.wx.imshow(
            np.abs(self.WaveletSignalPDV),
            extent=extent,
            cmap='PRGn',
            aspect='auto',
            origin='lower',  # pour que les basses fréquences soient en bas
            vmax=np.abs(self.WaveletSignalPDV).max()
            )


        self.wx.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        self.wx.set_ylabel("Frequency ("+str(self.Lbl_Freq_unit)+")")
        
        self.wx.set_title("Wavelet Spectrogram " + self.FName)

         # Restore previous zooming
        self.wx.set_xlim(wxlim)
        self.wx.set_ylim(wylim)

        # Save figure
        self.wfig.savefig(self.FName + '_SpectrogramWavelet.png')

        # Up date screen
        self.wcanvas.draw_idle()

    def ExtractWaveletVelocityNotebook(self):
        self.WVelocityProfile = []
        
        # Fclick on figure for value acquisition
        def onclick(event):
            if event.inaxes == self.wx:  # click on figure for value acquisition
                self.WVelocityProfile.append((event.xdata, event.ydata))
                print(f"Added point : {event.xdata:.4f}, {event.ydata:.4f}")
                self.wx.plot(event.xdata, event.ydata, 'rx')
                self.wcanvas.draw_idle()
    
        # Fonction de fin d'enregistrement : créer un nouvel onglet
        def stop_recording():
            self.wcanvas.mpl_disconnect(self.Wcid)
            self.Btn_ManualExtractWL.configure(text="Start extraction", command=self.ExtractWaveletVelocityNotebook)
            print("Extraction is over")
		
            # === Créer un nouvel onglet pour afficher les points ===
            self.Wframe_velocity = ttk.Frame(self.notebook)
            self.notebook.add(self.Wframe_velocity, text="Wavelet Velocity Extraction")
            self.notebook.select(self.Wframe_velocity)
    
            # Figure vide pour affichage des points extraits
            Wfig_vel, Wx_vel = plt.subplots(figsize=(3, 2))
            Wcanvas_vel = FigureCanvasTkAgg(Wfig_vel, master=self.Wframe_velocity)
            Wcanvas_vel.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(Wcanvas_vel, self.Wframe_velocity)
            toolbar.update()
            toolbar.pack(side=tk.TOP, fill=tk.X)
            for fig_num in plt.get_fignums():
                plt.close(fig_num)
    
            Wx_vel.set_title("Velocity profile " + self.FName)
            Wx_vel.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
            Wx_vel.set_ylabel("Velocity (m/s)")
            Wx_vel.grid(True)
    
            # Tracer les points extraits
            if self.WVelocityProfile:
                x, y = zip(*self.WVelocityProfile)
                x = np.asarray(list(x))
                y = np.asarray(list(y))
                # y2 = y*self.Freq_OM
                # y3 = y2 - self.Freq_pivot
                # y4 = y3*self.PDVFactor
                Wx_vel.plot(x, y, 'b.-')
                Wx_vel.legend()
                print ('Save velocity profile in '+self.FName+"VelocityProfileWavelet.csv")
                np.savetxt(self.FName+"VelocityProfileWavelet.csv", np.vstack((x ,y)).T, delimiter=',')
                print ('Save velocity plat in '+self.FName+"VelocityProfileWavelet.png")
                Wfig_vel.savefig(self.FName+'VelocityWavelet.png')
    
            Wcanvas_vel.draw_idle()
            
        self.Btn_ManualExtractWL.configure(text="Stop extraction", command=stop_recording)
        
        # Connexion du clic
        self.Wcid = self.wcanvas.mpl_connect('button_press_event', onclick)
        
    
    def AutoExtractVelocity_WL(self):
        # Load frequency boundaries
        FminWL = self.EntMinFreqWL.get()
        FmaxWL = self.EntMaxFreqWL.get()
        
        # Initiate error label
        self.Lbl_AutoExtract_ErrorLblWL.config(text = "")
        txtlblextract = self.Lbl_AutoExtract_ErrorLblWL.cget("text")
        
        if list(FminWL) == []:
            FminWL = np.min(self.WaveletFrequencies)/self.Freq_OM
            IndFminWL = np.argmin(np.abs(self.WaveletFrequencies - float(FminWL)*self.Freq_OM))
        else:
            if not(FminWL.replace('.','',1).isdigit()):
                txtlblextract = "Min frequency is not a number"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            IndFminWL = np.argmin(np.abs(self.WaveletFrequencies - float(FminWL)*self.Freq_OM))
        
        if list(FmaxWL) == []:
            FmaxWL = np.max(self.WaveletFrequencies)/self.Freq_OM
            IndFmaxWL = np.argmin(np.abs(self.WaveletFrequencies - float(FmaxWL)*self.Freq_OM))
        else:
            if not(FmaxWL.replace('.','',1).isdigit()):
                txtlblextract = "Max frequency is not a number"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            if float(FminWL)>float(FmaxWL):
                txtlblextract = "Min frequency is higher than max one"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            IndFmaxWL = np.argmin(np.abs(self.WaveletFrequencies - float(FmaxWL)*self.Freq_OM))
        
        # Load time boundaries
        TminWL = self.EntMinTimeWL.get()
        TmaxWL = self.EntMaxTimeWL.get()
        
        # Check min time boundary is a number, if yes, locate closest time in time vector
        if list(TminWL) == []:
            TminWL = np.min(self.Time)/self.Time_OM_Long
            IndTminWL = np.argmin(np.abs(self.Time - float(TminWL)*self.Time_OM_Long))
        else:
            if not(TminWL.replace('.','',1).isdigit()):
                txtlblextract = "Min time is not a number"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            IndTminWL = np.argmin(np.abs(self.Time - float(TminWL)*self.Time_OM_Long))
        
        # Check min time boundary is a number, if yes, locate closest time in time vector
        if list(TmaxWL) == []:
            TmaxWL = np.max(self.Time)/self.Time_OM_Long
            IndTmaxWL = np.argmin(np.abs(self.Time - float(TmaxWL)*self.Time_OM_Long))
        else:
            if not(TmaxWL.replace('.','',1).isdigit()):
                txtlblextract = "Max time is not a number"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            if float(TminWL)>float(TmaxWL):
                txtlblextract = "Min time is higher than max one"
                self.Lbl_AutoExtract_ErrorLblWL.config(text = txtlblextract)
                return
            IndTmaxWL = np.argmin(np.abs(self.Time - float(TmaxWL)*self.Time_OM_Long))
        
        self.WaveletSignalPDV_cut = self.WaveletSignalPDV[IndFmaxWL:IndFminWL+1, IndTminWL:(IndTmaxWL+1)]
        self.WaveletFrequencies_cut = self.WaveletFrequencies[IndFmaxWL:IndFminWL+1]
        self.Time_cut = self.Time[IndTminWL:(IndTmaxWL+1)]
        
        seuil = 0.12
        VecTime_cut = []
        Vec_IndProfFmaxWL = []
        Vec_Freq_IndProfFmaxWL = []
        F_BoundSup = []
        F_BoundInf = []
        for k in range(len(self.Time_cut)):
            if not(np.abs(np.max(self.WaveletSignalPDV_cut[:, k]))<seuil):
                VecTime_cut.append(self.Time_cut[k])
                Vec_IndProfFmaxWL.append(np.argmax(np.abs(self.WaveletSignalPDV_cut[:, k])))
                Vec_Freq_IndProfFmaxWL.append(self.WaveletFrequencies_cut[np.argmax(np.abs(self.WaveletSignalPDV_cut[:, k]))])
                
                ProfT_norm_tmp = np.abs(self.WaveletSignalPDV_cut[:, k])/np.max(np.abs(self.WaveletSignalPDV_cut[:, k]))
                ProfT_norm_tmp_Inf = ProfT_norm_tmp[(Vec_IndProfFmaxWL[-1]+1):]
                ProfT_norm_tmp_Sup = ProfT_norm_tmp[:Vec_IndProfFmaxWL[-1]]
                
                # Inferieur ou superieur a la limite fixee
                Bound_ProfWL = 0.5
                Sign_ProfT_norm_Inf_tmp = np.sign(ProfT_norm_tmp_Inf - Bound_ProfWL) + 1
                Sign_ProfT_norm_Sup_tmp = np.sign(ProfT_norm_tmp_Sup - Bound_ProfWL) + 1
                
                # Localisation de tous les non-zeros (points superieurs a la limite)
                Ind_ProfT_norm_Inf_tmp_SupBound = np.nonzero(Sign_ProfT_norm_Inf_tmp)
                Ind_ProfT_norm_Sup_tmp_SupBound = np.nonzero(Sign_ProfT_norm_Sup_tmp)
                
                # Intervalles de points encadrant le passage de la limite
                Inter_Ind_ProfT_norm_Inf_tmp_SupBound = ProfT_norm_tmp_Inf[Ind_ProfT_norm_Inf_tmp_SupBound[0][-1]:(Ind_ProfT_norm_Inf_tmp_SupBound[0][-1]+2)]
                Inter_Ind_ProfT_norm_Sup_tmp_SupBound = ProfT_norm_tmp_Sup[Ind_ProfT_norm_Sup_tmp_SupBound[0][0]-1:Ind_ProfT_norm_Sup_tmp_SupBound[0][0]+1]
                
                self.WaveletFrequencies_cutInf = self.WaveletFrequencies_cut[Vec_IndProfFmaxWL[-1]+1+Ind_ProfT_norm_Inf_tmp_SupBound[0][-1]:(Vec_IndProfFmaxWL[-1]+1+Ind_ProfT_norm_Inf_tmp_SupBound[0][-1]+2)]
                self.WaveletFrequencies_cutSup = self.WaveletFrequencies_cut[Ind_ProfT_norm_Sup_tmp_SupBound[0][0]-1:Ind_ProfT_norm_Sup_tmp_SupBound[0][0]+1]
                
                PosClosest_Ind_ProfT_norm_Inf_tmp = np.argmin(np.abs(Inter_Ind_ProfT_norm_Inf_tmp_SupBound-Bound_ProfWL))
                PosClosest_Ind_ProfT_norm_Sup_tmp = np.argmin(np.abs(Inter_Ind_ProfT_norm_Sup_tmp_SupBound-Bound_ProfWL))
                
                F_BoundSup.append(self.WaveletFrequencies_cutSup[PosClosest_Ind_ProfT_norm_Sup_tmp])
                F_BoundInf.append(self.WaveletFrequencies_cutInf[PosClosest_Ind_ProfT_norm_Inf_tmp])
            
        NameTabVelProfWL = "Velocity Profile WL"
        
        for tab_id in self.notebook.tabs():
            tab_text_tmp = self.notebook.tab(tab_id, "text")
            if tab_text_tmp == NameTabVelProfWL:
                self.notebook.forget(tab_id)
        
        self.Tab_VelProf_WL = ttk.Frame(self.notebook)
        self.notebook.add(self.Tab_VelProf_WL, text=NameTabVelProfWL)
        self.notebook.select(self.Tab_VelProf_WL)
        
        # Figure vide
        fig_VelProfWL, ax_VelProfWL = plt.subplots(figsize=(3, 2))
        canvas_VelProfWL = FigureCanvasTkAgg(fig_VelProfWL, master=self.Tab_VelProf_WL)
        canvas_VelProfWL.get_tk_widget().pack(side=tk.TOP, fill=None, expand=False)
        toolbarRWL = NavigationToolbar2Tk(canvas_VelProfWL, self.Tab_VelProf_WL)
        toolbarRWL.update()
        toolbarRWL.pack(side=tk.TOP, fill=tk.X)
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        
        ax_VelProfWL.set_title("Velocity profile WL " + self.FName)
        # ax_VelProfWL.set_xlabel("Time (s)")
        ax_VelProfWL.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        ax_VelProfWL.set_ylabel("Velocity (m/s)")
        ax_VelProfWL.grid(True)
        
        ax_VelProfWL.plot(np.asarray(VecTime_cut)/self.Time_OM_Long, np.asarray(Vec_Freq_IndProfFmaxWL)*self.PDVFactor, 'b.', label='Freq max velocity')
        ax_VelProfWL.fill_between(np.asarray(VecTime_cut)/self.Time_OM_Long, np.asarray(F_BoundInf)*self.PDVFactor, np.asarray(F_BoundSup)*self.PDVFactor, alpha=.3, linewidth=0, color='blue', label='+/-50% max velocity')
        ax_VelProfWL.legend()
        
        fig_VelProfWL.savefig(self.FName + '_AutoProfVelWL.png', dpi='figure')
        np.savetxt(self.FName+"VelocityProfile_WL.csv", np.vstack((VecTime_cut ,np.asarray(Vec_Freq_IndProfFmaxWL)*self.PDVFactor, np.asarray(F_BoundInf)*self.PDVFactor, np.asarray(F_BoundSup)*self.PDVFactor,)).T, delimiter=',')
        
        fig_VelProfWL.set_size_inches(3, 2)
        canvas_VelProfWL.draw_idle()
        
    
    def CreateSTFTPDVInteractive(self, parent):
        #Give value for initial slider
        self.param = self.nperseg
    
        self.fig, self.ax = plt.subplots(figsize=(3, 2))
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
    		
        # self.label1 = tk.Label(parent, text=f"STFT Window (nperseg) : {self.nperseg} pt , {self.WindowsSize*1e9:.3f} ns, Window : {self.STFTPDVWindow}")
        self.label1 = tk.Label(parent, text=f"STFT Window (nperseg) : {self.nperseg} pt ,"+ str(self.WindowsSize*self.Freq_OM) + self.Lbl_Time_unit_Court + ", Window : {self.STFTPDVWindow}")
        self.label1.pack()
        
        # self.label2 = tk.Label(parent, text=f"Number of points : {len(self.Time)} pt, FAcquisition (GS/s): {self.FAcquisition*1e-9:e}")
        self.label2 = tk.Label(parent, text=f"Number of points : {len(self.Time)} pt, FAcquisition ("+str(self.Lbl_Freq_unit)+"): "+str(self.FAcquisition/self.Freq_OM))
        self.label2.pack()
            
        self.ax.set_title("Spectrogram + self.FName")
        # self.ax.set_xlabel("Time (s)") #ICITIME
        # self.ax.set_ylabel("Frequency (Hz)") #ICITIME
        self.ax.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        self.ax.set_ylabel("Frequency ("+str(self.Lbl_Freq_unit)+")")
        
        self.STFTPDVWindow_var = tk.StringVar(value=self.STFTPDVWindow)
        fenetres = ['hann', 'hamming', 'blackman', 'bartlett', 'flattop']

        self.LblTitle_ParamSTFT = tk.Label(parent, text="STFT Parameters", font=("Arial", 14, "bold"))
        self.LblTitle_ParamSTFT.pack(anchor = "w")
        
        self.Frame_WindowTimeBase = ttk.Frame(parent)
        self.Frame_WindowTimeBase.pack(anchor="w")
        
        self.Frame_ParamSTFT = ttk.Frame(self.Frame_WindowTimeBase)
        self.Frame_ParamSTFT.pack(anchor="w")
        
        ttk.Label(self.Frame_ParamSTFT, text="Windows STFT :").pack(anchor = "w", padx=5, side=tk.LEFT)
        
        combo = ttk.Combobox(self.Frame_ParamSTFT, textvariable=self.STFTPDVWindow_var, values=fenetres, state="readonly")
        combo.pack(anchor = "w", padx=5, side=tk.LEFT)
        combo.bind('<<ComboboxSelected>>', lambda e: self.update_STFTPDVInteractiveplot(self.slider.get()))
        
        self.Btnnperseg = tk.Button(self.Frame_ParamSTFT, text="nperseg", command=self.Updatenperseg)
        self.Btnnperseg.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.Entnperseg = ttk.Entry(self.Frame_ParamSTFT, width = 5)
        self.Entnperseg.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.slider = ttk.Scale(self.Frame_ParamSTFT, from_=2, to=2048, orient='horizontal')
        self.slider.set(self.nperseg)
        self.slider.pack(anchor = "w", padx=5, side=tk.LEFT)
        self.slider.configure(command=self.update_STFTPDVInteractiveplot)
        
        self.BtnHelp_Baseline = tk.Button(self.Frame_ParamSTFT, text="?", command=self.InterfaceHelpBaseline)
        self.BtnHelp_Baseline.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        ttk.Label(self.Frame_ParamSTFT, text="Baseline management :").pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.BaseLineManag = tk.Button(self.Frame_ParamSTFT, text="Delete", command = self.BaseLineDelete)
        self.BaseLineManag.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        ###
        self.Frame_ManualExtractSTFT = ttk.Frame(parent)
        self.Frame_ManualExtractSTFT.pack(anchor="w")
        
        self.BtnHelp_ManualExtractSTFT = tk.Button(self.Frame_ManualExtractSTFT, text="?", command=self.InterfaceHelpManVelExtr)
        self.BtnHelp_ManualExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.LblTitle_ManualExtractSTFT = tk.Label(self.Frame_ManualExtractSTFT, text="Manual Velocity Extraction", font=("Arial", 14, "bold"))
        self.LblTitle_ManualExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.Frame_ManualExtractSTFT2 = ttk.Frame(parent)
        self.Frame_ManualExtractSTFT2.pack(anchor="w")
        
        self.Instruction_ManualExtractSTFT = tk.Label(self.Frame_ManualExtractSTFT2, text="Manually extract velocity by clicking on the spectrogram.")
        self.Instruction_ManualExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.Velocity_buttonSTFT = tk.Button(self.Frame_ManualExtractSTFT2, text="Start Extraction", command=self.ExtractVelocityNotebook)
        self.Velocity_buttonSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        ###
        self.Frame_AutoExtractSTFT = ttk.Frame(parent)
        self.Frame_AutoExtractSTFT.pack(anchor="w", pady=5)
        
        self.BtnHelp_AutoExtractSTFT = tk.Button(self.Frame_AutoExtractSTFT, text="?", command=self.InterfaceHelpAutoVelExtr)
        self.BtnHelp_AutoExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.LblTitle_AutoExtractSTFT = tk.Label(self.Frame_AutoExtractSTFT, text="Automatic Velocity Extraction", font=("Arial", 14, "bold"))
        self.LblTitle_AutoExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.Frame_AutoExtractSTFT2 = ttk.Frame(parent)
        self.Frame_AutoExtractSTFT2.pack(anchor="w", pady=5)
        
        self.Instruction_AutoExtractSTFT = tk.Label(self.Frame_AutoExtractSTFT2, text="Automatic extraction in a defined window")
        self.Instruction_AutoExtractSTFT.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        MinFreqAutoVelExtract_frame = ttk.Frame(self.Frame_AutoExtractSTFT2)
        MinFreqAutoVelExtract_frame.pack(anchor = "w", padx=5, side=tk.LEFT)
        # ttk.Label(MinFreqAutoVelExtract_frame, text="Min Freq (GHz)").pack(anchor="w") #ICITIME
        ttk.Label(MinFreqAutoVelExtract_frame, text="Min Freq ("+str(self.Lbl_Freq_unit)+")").pack(anchor="w")
        self.EntMinFreq = tk.Entry(MinFreqAutoVelExtract_frame, width=15)
        self.EntMinFreq.pack()
        
        MaxFreqAutoVelExtract_frame = ttk.Frame(self.Frame_AutoExtractSTFT2)
        MaxFreqAutoVelExtract_frame.pack(anchor = "w", padx=5, side=tk.LEFT)
        # ttk.Label(MaxFreqAutoVelExtract_frame, text="Max Freq (GHz)").pack(anchor="w") #ICITIME
        ttk.Label(MaxFreqAutoVelExtract_frame, text="Max Freq ("+str(self.Lbl_Freq_unit)+")").pack(anchor="w")
        self.EntMaxFreq = tk.Entry(MaxFreqAutoVelExtract_frame, width=15)
        self.EntMaxFreq.pack()
        
        MinTimeAutoVelExtract_frame = ttk.Frame(self.Frame_AutoExtractSTFT2)
        MinTimeAutoVelExtract_frame.pack(anchor = "w", padx=5, side=tk.LEFT)
        # ttk.Label(MinTimeAutoVelExtract_frame, text="T min (µs)").pack(anchor="w") #ICITIME
        ttk.Label(MinTimeAutoVelExtract_frame, text="T min ("+str(self.Lbl_Time_unit_Long)+")").pack(anchor="w")
        self.EntMinTime = tk.Entry(MinTimeAutoVelExtract_frame, width=15)
        self.EntMinTime.pack()
        
        MaxTimeAutoVelExtract_frame = ttk.Frame(self.Frame_AutoExtractSTFT2)
        MaxTimeAutoVelExtract_frame.pack(anchor = "w", padx=5, side=tk.LEFT)
        # ttk.Label(MaxTimeAutoVelExtract_frame, text="T max (µs)").pack(anchor="w") #ICITIME
        ttk.Label(MaxTimeAutoVelExtract_frame, text="T max ("+str(self.Lbl_Time_unit_Long)+")").pack(anchor="w")
        self.EntMaxTime = tk.Entry(MaxTimeAutoVelExtract_frame, width=15)
        self.EntMaxTime.pack()
        
        AutoVelocity_button = tk.Button(self.Frame_AutoExtractSTFT2, text="Automatic extraction", command=self.ExtractVelocityNotebookAuto)
        AutoVelocity_button.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        self.ErrorAutoVel_Lbl = ttk.Label(self.Frame_AutoExtractSTFT2, text="", foreground="red")
        self.ErrorAutoVel_Lbl.pack(anchor = "w", padx=5, side=tk.LEFT)
        
        # Calcul initial
        self.SetSTFTPDV(self.nperseg)
        
        # Initial Plot
        self.quadmesh = self.ax.pcolormesh(
            self.Time_stft/self.Time_OM_Long,
            self.FePDV/self.Freq_OM,
            np.abs(self.PDVSpectrogram),
            shading='gouraud'
            )
        
        # Affichage initial
        self.PDVSpectrogramActive = np.abs(self.PDVSpectrogram)
        # self.quadmesh = self.ax.pcolormesh(self.Time_stft, self.FePDV, self.PDVSpectrogramActive, shading='gouraud') #ICITIME
        self.quadmesh = self.ax.pcolormesh(self.Time_stft/self.Time_OM_Long, self.FePDV/self.Freq_OM, self.PDVSpectrogramActive, shading='gouraud')
        
        self.ax.set_ylim(min(self.FePDV/self.Freq_OM), max(self.FePDV/self.Freq_OM))
        self.ax.set_title("Spectrogram " + self.FName)
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        self.canvas.draw_idle()
     

    def BaseLineDelete(self):
        self.BaseLineManag.configure(text="Reset")
        self.BaseLineManag.configure(command=self.ResetBaseline)
        BaseLineFreq = np.argmax(self.PDVSpectrogram[:, 5])
        VecBaseLine = np.abs(self.PDVSpectrogram[:, 5])
        for k in range(len(self.Time_stft)):
            self.PDVSpectrogramActive[:, k] = np.abs(self.PDVSpectrogram[:, k]) - VecBaseLine*np.abs(self.PDVSpectrogram[BaseLineFreq, k])/np.abs(self.PDVSpectrogram[BaseLineFreq, 5])
        self.quadmesh = self.ax.pcolormesh(self.Time_stft, self.FePDV, self.PDVSpectrogramActive, shading='gouraud')
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        self.fig.tight_layout()
        self.canvas.draw_idle()
        self.fig.savefig(self.FName+'Spectrogram.png', dpi=200)
        
    def ResetBaseline(self):
        self.BaseLineManag.configure(text="Delete")
        self.BaseLineManag.configure(command = self.BaseLineDelete)
        self.PDVSpectrogramActive = np.abs(self.PDVSpectrogram)
        self.quadmesh = self.ax.pcolormesh(self.Time_stft/self.Time_OM_Long, self.FePDV/self.Freq_OM, self.PDVSpectrogramActive, shading='gouraud')
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        self.canvas.draw_idle()
        
    
    def InterfaceHelpManVelExtr(self):
        self.Inter_Help_ManVelExtract = tk.Tk()
        
        self.LblTitle_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Manual Extraction Velocity", font=("Arial", 14, "bold"))
        self.LblTitle_Help_ManVelExtract.pack(anchor="w", pady=5)
        
        self.Lbl1_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Once clicked on 'Start Extraction', you can click on the spectrogram.")
        self.Lbl1_Help_ManVelExtract.pack(anchor="w")
        
        self.Lbl2_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Clicked position (time and frequency) will be saved and a red cross added to the spectrogram.")
        self.Lbl2_Help_ManVelExtract.pack(anchor="w")
        
        self.Lbl3_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Click on 'Stop Extraction' to end collecting points.")
        self.Lbl3_Help_ManVelExtract.pack(anchor="w")
        
        self.Lbl4_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="All the points will be displayed in their own figure in a new tab.")
        self.Lbl4_Help_ManVelExtract.pack(anchor="w")
        
        self.Lbl5_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Data points are also save in a csv file labelled '{Initital Shot File Name} + VelocityProfile.csv'.")
        self.Lbl5_Help_ManVelExtract.pack(anchor="w")
        
        self.Lbl6_Help_ManVelExtract = tk.Label(self.Inter_Help_ManVelExtract, text="Clicking on 'Stop Extraction' erase previous, if any, figure and csv file.")
        self.Lbl6_Help_ManVelExtract.pack(anchor="w")
        
        self.Btn_QuitHelp_ManVelExtract = tk.Button(self.Inter_Help_ManVelExtract, text="Leave", command=lambda: self.Inter_Help_ManVelExtract.destroy())
        self.Btn_QuitHelp_ManVelExtract.pack(anchor="w")
        
        self.Inter_Help_ManVelExtract.mainloop()
        
    
    def InterfaceHelpAutoVelExtr(self):
        self.Inter_Help_AutoVelExtract = tk.Tk()
        
        self.LblTitle_Help_AutoVelExtract = tk.Label(self.Inter_Help_AutoVelExtract, text="Automatic Extraction Velocity", font=("Arial", 14, "bold"))
        self.LblTitle_Help_AutoVelExtract.pack(anchor="w", pady=5)
        
        self.Lbl1_Help_AutoVelExtract = tk.Label(self.Inter_Help_AutoVelExtract, text="For each time step, extract frequency corresponding to the maximal amplitude.")
        self.Lbl1_Help_AutoVelExtract.pack(anchor="w")
        
        self.Lbl2_Help_AutoVelExtract = tk.Label(self.Inter_Help_AutoVelExtract, text="Frequency interval has to be given.")
        self.Lbl2_Help_AutoVelExtract.pack(anchor="w")
        
        self.Lbl3_Help_AutoVelExtract = tk.Label(self.Inter_Help_AutoVelExtract, text="Time window can be given. If not, extraction is done for the entire spectrogram.")
        self.Lbl3_Help_AutoVelExtract.pack(anchor="w")
        
        self.Btn_QuitHelp_AutoVelExtract = tk.Button(self.Inter_Help_AutoVelExtract, text="Leave", command=lambda: self.Inter_Help_AutoVelExtract.destroy())
        self.Btn_QuitHelp_AutoVelExtract.pack(anchor="w")
        
        self.Inter_Help_AutoVelExtract.mainloop()
        
    
    def InterfaceHelpBaseline(self):
        self.Inter_Help_Baseline = tk.Tk()
        
        self.LblTitle_Help_Baseline = tk.Label(self.Inter_Help_Baseline, text="Baseline management", font=("Arial", 14, "bold"))
        self.LblTitle_Help_Baseline.pack(anchor="w", pady=5)
        
        self.Lbl1_Help_Baseline = tk.Label(self.Inter_Help_Baseline, text="Extract a reference frequency spectrum at a time step before first wave arrival.")
        self.Lbl1_Help_Baseline.pack(anchor="w")
        
        self.Lbl2_Help_Baseline = tk.Label(self.Inter_Help_Baseline, text="Define amplitude at pivot frequency.")
        self.Lbl2_Help_Baseline.pack(anchor="w")
        
        self.Lbl3_Help_Baseline = tk.Label(self.Inter_Help_Baseline, text="Reference frequency spectrum is substracted at each time step once normalized by amplitude at pivot frequency for this time step.")
        self.Lbl3_Help_Baseline.pack(anchor="w")
        
        self.Lbl4_Help_Baseline = tk.Label(self.Inter_Help_Baseline, text="Reseting baseline reload initially calculated spectrogram.")
        self.Lbl4_Help_Baseline.pack(anchor="w")
        
        self.Btn_QuitHelp_Baseline = tk.Button(self.Inter_Help_Baseline, text="Leave", command=lambda: self.Inter_Help_Baseline.destroy())
        self.Btn_QuitHelp_Baseline.pack(anchor="w")
        
        self.Inter_Help_Baseline.mainloop()
        
    
    def Updatenperseg(self):
		try:
        	# Read width
            nperseg_str = self.nperseg_var.get()
            self.nperseg = int(float(nperseg_str))
            print ("self.nperseg :",self.nperseg )
            if self.nperseg <= 2:
                print("nperseg need to be >2.")
                return
        except ValueError:
            print("Invalid value")
            return
					
        Newnperseg = float(self.Entnperseg.get())
        self.update_STFTPDVInteractiveplot(Newnperseg)
        
    
    def update_STFTPDVInteractiveplot(self, val):
        
        val = float(val)
        window = self.STFTPDVWindow_var.get()
        self.STFTPDVWindow = window
        self.nperseg=int(val)
        
        self.label1.config(
            text=f"STFT Window (nperseg) : {self.nperseg} pt ,"+ str(self.WindowsSize*self.Freq_OM) + self.Lbl_Time_unit_Court + f", Window : {self.STFTPDVWindow}")
        
        self.label2.config(
            text=f"Number of points : {len(self.Time)} pt, FAcquisition ("+str(self.Lbl_Freq_unit)+"): "+str(self.FAcquisition/self.Freq_OM))
        
        # save zooming
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        # Update STFT
        self.SetSTFTPDV(self.nperseg)

        # Recalculer STFT
        self.SetSTFTPDV(self.param)
        self.PDVSpectrogramActive = np.abs(self.PDVSpectrogram)
        
        if self.BaseLineManag.cget('text')=="Reset":
            self.BaseLineDelete()
    
        # Effacer seulement le contenu des axes
        self.ax.clear()
        
        # Initial Plot
        self.quadmesh = self.ax.pcolormesh(
            self.Time_stft/self.Time_OM_Long,
            self.FePDV/self.Freq_OM,
            np.abs(self.PDVSpectrogram),
            shading='gouraud'
            )
        
        # self.ax.set_xlabel("Time (s)")
        # self.ax.set_ylabel("Frequency (Hz)")
        self.ax.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        self.ax.set_ylabel("Frequency ("+str(self.Lbl_Freq_unit)+")")
        self.ax.set_title("Spectrogram " + self.FName)
        #go to previous zooming
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        #save figures
        self.fig.savefig(self.FName+'Spectrogram.png')
        # upate figures
        self.canvas.draw_idle()
    
    # Fonction d'extraction manuelle du profil de vitesse
    def ExtractVelocityNotebook(self):
        self.VelocityProfile = []
        
        # Fonction de clic dans le spectrogramme interactif
        def onclick(event):
            if event.inaxes == self.ax:
                # self.VelocityProfile.append((event.xdata, event.ydata*self.PDVFactor))
                self.VelocityProfile.append((event.xdata, event.ydata))
                print(f"addet point : {event.xdata:.4f}, {event.ydata:.4f}")
                self.ax.plot(event.xdata, event.ydata, 'rx')
                self.canvas.draw_idle()
        
        # Fonction de fin d'enregistrement : créer un nouvel onglet
        def stop_recording():
            self.canvas.mpl_disconnect(self.cid)
            print("Extraction is over")
            
            self.Velocity_buttonSTFT.configure(text="Start extraction", command=self.ExtractVelocityNotebook)
            
            # === Créer un nouvel onglet pour afficher les points ===

            self.frame_velocity = ttk.Frame(self.notebook)
            self.notebook.add(self.frame_velocity, text="STFT Velocity Extraction")
            self.notebook.select(self.frame_velocity)
            
            self.Ind_Freq_pivot = np.argmax(self.PDVSpectrogram[:, 5])
            self.Freq_pivot = self.FePDV[self.Ind_Freq_pivot]
            
            fig_vel, ax_vel = plt.subplots(figsize=(3, 2))
            canvas_vel = FigureCanvasTkAgg(fig_vel, master=self.frame_velocity)
            canvas_vel.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas_vel, self.frame_velocity)
            toolbar.update()
            toolbar.pack(side=tk.TOP, fill=tk.X)
            for fig_num in plt.get_fignums():
                plt.close(fig_num)
            
            ax_vel.set_title("Velocity profile " + self.FName)
            ax_vel.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
            ax_vel.set_ylabel("Velocity (m/s)") #ICITIME
            
            ax_vel.grid(True)
            if self.VelocityProfile:
                x, y = zip(*self.VelocityProfile)
                x = np.asarray(list(x))
                y = np.asarray(list(y))
                y2 = y*self.Freq_OM
                y3 = y2 - self.Freq_pivot
                y4 = y3*self.PDVFactor
                ax_vel.plot(x, y4, 'rx-')
                ax_vel.legend()
                print ('Save velocity profile in '+self.FName+"VelocityProfile.csv")
                np.savetxt(self.FName+"VelocityProfile.csv", np.vstack((x ,y)).T, delimiter=',')
                print ('Save velocity plat in '+self.FName+"VelocityProfile.png")
                fig_vel.savefig(self.FName+'Velocity.png')
        self.Velocity_buttonSTFT.configure(text="Stop extraction", command=stop_recording)
        
        # Connexion click
        self.cid = self.canvas.mpl_connect('button_press_event', onclick)
        
    
    # Fonction d'extraction automatique du profil de vitesse
    def ExtractVelocityNotebookAuto(self):
        Fmin = self.EntMinFreq.get()
        Fmax = self.EntMaxFreq.get()
        
        self.ErrorAutoVel_Lbl.config(text = "")
        txtlblextract = self.ErrorAutoVel_Lbl.cget("text")
        
        if list(Fmin) == []:
            # Fmin = np.min(self.FePDV)*1e-9 #ICITIME
            Fmin = np.min(self.FePDV)/self.Freq_OM
            # IndFmin = np.argmin(np.abs(self.FePDV - float(Fmin)*1e9)) #ICITIME
            IndFmin = np.argmin(np.abs(self.FePDV - float(Fmin)*self.Freq_OM))
        else:
            if not(Fmin.replace('.','',1).isdigit()):
                txtlblextract = "Min frequency is not a number"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            # IndFmin = np.argmin(np.abs(self.FePDV - float(Fmin)*1e9)) #ICITIME
            IndFmin = np.argmin(np.abs(self.FePDV - float(Fmin)*self.Freq_OM))
        
        if list(Fmax) == []:
            # Fmax = np.max(self.FePDV)*1e-9
            # IndFmax = np.argmin(np.abs(self.FePDV - float(Fmax)*1e9)) #ICITIME
            Fmax = np.max(self.FePDV)/self.Freq_OM
            IndFmax = np.argmin(np.abs(self.FePDV - float(Fmax)*self.Freq_OM))
        else:
            if not(Fmax.replace('.','',1).isdigit()):
                txtlblextract = "Max frequency is not a number"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            if float(Fmin)>float(Fmax):
                txtlblextract = "Min frequency is higher than max one"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            # IndFmax = np.argmin(np.abs(self.FePDV - float(Fmax)*1e9)) #ICITIME
            IndFmax = np.argmin(np.abs(self.FePDV - float(Fmax)*self.Freq_OM))
        
        Tmin = self.EntMinTime.get()
        Tmax = self.EntMaxTime.get()
        
        # Check min time boundary is a number, if yes, locate closest time in time vector
        if list(Tmin) == []:
            # Tmin = np.min(self.Time_stft)*1e6 #ICITIME
            Tmin = np.min(self.Time_stft)/self.Time_OM_Long
            # IndTmin = np.argmin(np.abs(self.Time_stft - float(Tmin)*1e-6)) #ICITIME
            IndTmin = np.argmin(np.abs(self.Time_stft - float(Tmin)*self.Time_OM_Long))
        else:
            if not(Tmin.replace('.','',1).isdigit()):
                txtlblextract = "Min time is not a number"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            # IndTmin = np.argmin(np.abs(self.Time_stft - float(Tmin)*1e-6)) #ICITIME
            IndTmin = np.argmin(np.abs(self.Time_stft - float(Tmin)*self.Time_OM_Long))
        
        # Check min time boundary is a number, if yes, locate closest time in time vector
        if list(Tmax) == []:
            # Tmax = np.max(self.Time_stft)*1e6 #ICITIME
            # IndTmax = np.argmin(np.abs(self.Time_stft - float(Tmax)*1e-6)) #ICITIME
            Tmax = np.max(self.Time_stft)/self.Time_OM_Long
            IndTmax = np.argmin(np.abs(self.Time_stft - float(Tmax)*self.Time_OM_Long))
        else:
            if not(Tmax.replace('.','',1).isdigit()):
                txtlblextract = "Max time is not a number"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            if float(Tmin)>float(Tmax):
                txtlblextract = "Min time is higher than max one"
                self.ErrorAutoVel_Lbl.config(text = txtlblextract)
                return
            # IndTmax = np.argmin(np.abs(self.Time_stft - float(Tmax)*1e-6)) #ICITIME
            IndTmax = np.argmin(np.abs(self.Time_stft - float(Tmax)*self.Time_OM_Long))
        
        self.PDVSpectrogram_cut = self.PDVSpectrogramActive[IndFmin:IndFmax, IndTmin:IndTmax]
        self.FePDV_cut = self.FePDV[IndFmin:IndFmax]
        self.Time_stft_cut = self.Time_stft[IndTmin:IndTmax]
        
        IndFmaxT = np.argmax(self.PDVSpectrogram_cut, axis=0)
        self.Prof_FMax = self.FePDV_cut[IndFmaxT]
        
        self.Freq_pivot = self.Prof_FMax[0]
                
        Bound_Prof = 0.5
        Bound_v_Inf = []
        Bound_v_Sup = []
        
        for k in range((IndTmax-IndTmin)):
            Mat_PDVSpec_Cut_Norm = np.abs(self.PDVSpectrogram_cut[:, k])/np.max(np.abs(self.PDVSpectrogram_cut[:, k]))
            
            Mat_PDVSpec_Cut_Inf = Mat_PDVSpec_Cut_Norm[:IndFmaxT[k]]
            Mat_PDVSpec_Cut_Sup = Mat_PDVSpec_Cut_Norm[IndFmaxT[k]+1:]
            
            self.FePDV_Inf = self.FePDV_cut[:IndFmaxT[k]]
            self.FePDV_Sup = self.FePDV_cut[IndFmaxT[k]+1:]
            
            Sign_PDVSpec_Cut_Inf = np.sign(Mat_PDVSpec_Cut_Inf - Bound_Prof)
            Sign_PDVSpec_Cut_Sup = np.sign(Mat_PDVSpec_Cut_Sup - Bound_Prof)
            
            Zero_PDVSpec_Cut_Inf = np.nonzero(Sign_PDVSpec_Cut_Inf - 1)
            Zero_PDVSpec_Cut_Sup = np.nonzero(Sign_PDVSpec_Cut_Sup - 1)
            
            if np.sum(Sign_PDVSpec_Cut_Inf + 1) == 0:
                Bound_v_Inf.append(self.FePDV_Inf[-1])
            else:
                IndLstZero_Inf = Zero_PDVSpec_Cut_Inf[0][-1]
                if IndLstZero_Inf == len(Mat_PDVSpec_Cut_Inf):
                    Bound_v_Inf.append(self.FePDV_Inf[-1])
                else:
                    CloseBound_Inf = np.argmin(np.abs(Mat_PDVSpec_Cut_Inf[IndLstZero_Inf:IndLstZero_Inf+2] - Bound_Prof))
                    Bound_v_Inf.append(self.FePDV_Inf[IndLstZero_Inf+CloseBound_Inf])
            
            if np.sum(Sign_PDVSpec_Cut_Sup + 1) == 0:
                Bound_v_Sup.append(self.FePDV_Sup[0])
            else:
                IndFstZero_Inf = Zero_PDVSpec_Cut_Sup[0][0]
                if IndFstZero_Inf == 0 :
                    Bound_v_Sup.append(self.FePDV_Sup[0])
                else:
                    CloseBound_Sup = np.argmin(np.abs(Mat_PDVSpec_Cut_Sup[IndFstZero_Inf-1:IndFstZero_Inf+1] - Bound_Prof))
                    Bound_v_Sup.append(self.FePDV_Sup[IndFstZero_Inf + CloseBound_Sup])
        
        self.Prof_FMax = self.Prof_FMax - self.Freq_pivot
        Bound_v_Sup = Bound_v_Sup - self.Freq_pivot
        Bound_v_Inf = Bound_v_Inf - self.Freq_pivot
        
        self.Prof_FMax = self.Prof_FMax*self.PDVFactor
        Bound_v_Sup = np.asarray(Bound_v_Sup)*self.PDVFactor
        Bound_v_Inf = np.asarray(Bound_v_Inf)*self.PDVFactor
        
        NameOngProf = "Freq Profile"
        
        for tab_id in self.notebook.tabs():
            tab_text_tmp = self.notebook.tab(tab_id, "text")
            if tab_text_tmp == NameOngProf:
                self.notebook.forget(tab_id)
        
        self.frame_FreqProf = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_FreqProf, text=NameOngProf)
        self.notebook.select(self.frame_FreqProf)
        
        fig_velR, ax_velR = plt.subplots(figsize=(6, 4))
        canvas_velR = FigureCanvasTkAgg(fig_velR, master=self.frame_FreqProf)
        canvas_velR.get_tk_widget().pack(side=tk.TOP, fill=None, expand=False)
        toolbarR = NavigationToolbar2Tk(canvas_velR, self.frame_FreqProf)
        toolbarR.update()
        toolbarR.pack(side=tk.TOP, fill=tk.X)
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        
        ax_velR.set_title("Velocity profile " + self.FName)
        # ax_velR.set_xlabel("Time (s)")
        # ax_velR.set_ylabel("Velocity (m/s)") #ICITIME
        # ax_velR.set_xlabel("Time (s)")
        # ax_velR.set_ylabel("Velocity (m/s)")
        ax_velR.set_xlabel("Time ("+str(self.Lbl_Time_unit_Long)+")")
        ax_velR.set_ylabel("Velocity (m/"+str(self.Lbl_Time_unit_Long)+")")
        ax_velR.grid(True)
        
        ax_velR.plot(self.Time_stft[IndTmin:IndTmax]/self.Time_OM_Long, self.Prof_FMax, 'r.-', label='Max')
        ax_velR.fill_between(self.Time_stft[IndTmin:IndTmax]/self.Time_OM_Long, Bound_v_Inf, Bound_v_Sup, alpha=.3, linewidth=0, color='red', label='+/-50% max velocity')
        ax_velR.legend()
        
        fig_velR.savefig(self.FName + '_AutoProfVel.png', dpi='figure')
        
        fig_velR.set_size_inches(3, 2)
        canvas_velR.draw_idle()
        
        FileNameSave = "Save_AutoExtract_Velocity_Prof.csv"
        np.savetxt(FileNameSave, (self.Time_stft[IndTmin:IndTmax]/self.Time_OM_Long, self.Prof_FMax*self.PDVFactor, (np.asarray(Bound_v_Inf)*self.PDVFactor), (np.asarray(Bound_v_Sup)*self.PDVFactor)), header="time (s), max vel (m/s), max vel +3dB (m/s), max vel -3dB (m/s)", delimiter=',', newline=';')
        
    
    def NotebookGraphSpectrogram(self, parent):
        #raw datas plot
        fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(3, 2))
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2Tk(canvas, parent)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
        
        # axs[0].plot(self.Time, self.Tension) #ICITIME
        axs[0].plot(self.Time/self.Time_OM_Long, self.Tension)
        axs[0].set_title("Signal (t) " + self.FName)
        # axs[0].set_xlabel("Time (s)") #ICITIME
        axs[0].set_xlabel("Time ("+ str(self.Lbl_Time_unit_Long) + ")")
        axs[0].set_ylabel("Amplitude (V)")
        axs[0].grid()
        
        # axs[1].plot(self.PDVSignalFFTTime, np.abs(self.HSignalFFT)) #ICITIME
        axs[1].plot(self.PDVSignalFFTTime/self.Freq_OM, np.abs(self.HSignalFFT))
        axs[1].set_title("FFT")
        # axs[1].set_xlabel("Fe(Hz)")
        axs[1].set_xlabel("Fe ("+self.Lbl_Freq_unit+")")
        axs[1].set_ylabel("Magnitude")
        axs[1].set_yscale('log')
        # axs[1].set_xlim(0, self.ChainResponse*1.5)
        axs[1].set_xlim(0, self.ChainResponse*1.5/self.Freq_OM) #ICITIME
        axs[1].grid()
        
        fig.tight_layout()
        print ('Save Raw Data'+self.FName+'RawData.png')
        fig.savefig(self.FName+'RawData.png')
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        canvas.draw_idle()
        
    
    def runSTFTPDVInteractive(self):
        # Check if the OS is Windows
        if os.name == 'nt':
            self.root.state("zoomed")
        else:
                # For Unix-like systems, use `_NET_WM_STATE_MAXIMIZED_VERT` and `_NET_WM_STATE_MAXIMIZED_HORZ`
                # and configure the window to be maximized.
           self.root.attributes('-zoomed', True)
        self.root.mainloop()
      
    def PDVParameters(self):
        #calculation of pdv parameters on tab data & operation. 
        try:
            freq_ghz = float(self.ChainResponse_var.get())
            wavelength_nm = float(self.LambdaLaser_var.get())
            Shift_nm=float(self.Shift_var.get())
        except ValueError:
            print("Value error - check")
            return
        
        # SI Unit : GHz -> Hz, nm -> m
        freq_hz = freq_ghz
        wavelength_m = wavelength_nm
        Shift_nm=Shift_nm
        
        # Velocity calculation : v = f × λ / 2
        max_velocity = freq_hz * wavelength_nm / 2
        VPivot=Shift_nm*wavelength_nm/2
        self.MaxVelocityForChainResponse=max_velocity
        self.VPivot=VPivot
        self.MaxVelocityForChainResponse_var.set(f"{max_velocity:.2f}")
        self.VPivot_var.set(f"{self.VPivot:.2f}")
        
    
    def PDVReport(self): #Pdf report of shot
        
        print ("ReportVH pdf "+self.ShotNumber+' '+self. FName)
        w, h = A4
        c = canvas.Canvas("VHReport"+self.ShotNumber+".pdf", pagesize=A4)
        c.setFont("Helvetica", 8)
        c.drawString(10, h - 25, "Parameters : "+self.FName)
        c.drawString(10, h - 37, "Chain Response (GHz) : "+str(self.ChainResponse*1e-9) + "  >>>Max. Corresponding Velocity (m/s) :  "+str(self.ChainResponse*self.PDVFactor)) 
        c.drawString(10, h - 49, "PDV Shift (GHz)      : "+str(self.Shift*1e-9) + "  >>>Max. Corresponding Velocity (m/s)  :  "+str(self.Shift*self.PDVFactor))             
        c.drawString(10, h - 61, "FAquisition (GS/s)   :"+str(self.FAcquisition*1e-9)) 
        c.drawString(10, h - 73, "PDVFactor m/s/Hz :"+ str(self.PDVFactor)+"    Windows size (points) : "+str(self.SegSize)) 
        #c.drawString(10, h - 85, "BaseLine Freq_min (GHZ) :"\
        #             + str(self.FreqBaseLine_min*1e-9)\
        #             + " -  BaseLine Freq_max (GHZ) : "+str(self.FreqBaseLine_max*1e-9))

        img = ImageReader('PDVSpectrogram.png')
        c.drawImage(img, 10,h - 380, width=280,height=280)
        img2 = ImageReader("RawData.png")
        c.drawImage(img2, 300,h - 380, width=280,height=280)
    
        #img3 = ImageReader("PDVSignalSelected.png")
        #c.drawImage(img3, 10,h - 750, width=280,height=364)

        #img4 = ImageReader("PDVSignalWithoutBaseline.png")
        #c.drawImage(img4, 300,h - 750, width=280,height=364)

        c.showPage()
        c.save() 
    
        return
##########################TOOLS    
def PDVDesign(LambdaLaser,PDVShift,ChainResponse,TargetVelocity) :
    
    PDVShift=PDVShift*1e+9
    ChainResponse=ChainResponse*1e9
    print ("Chain Parameters **************************")
    print("Shift : (Hz) :",f"{PDVShift:e}")
    print("ChainResponse: (Hz) :",f"{ChainResponse:e}")
    print("TargetVelocity: (m/s) :",f"{TargetVelocity:e}")

    print ("Design*****")
    VPivot=LambdaLaser*PDVShift/2
    print("VPivot (m/s) :",VPivot)
    MaxVelocityForChainResponse=LambdaLaser*ChainResponse/2
    print("MaxVelocityForChainResponse (VP=0) (m/s) :",MaxVelocityForChainResponse)
    TargetVelocityFre=2/LambdaLaser*TargetVelocity
    TargetVelocityFreShift=2/LambdaLaser*abs(TargetVelocity-VPivot)
    print("TargetVelocity corres. Freq No shift: (GHz) :",f"{TargetVelocityFre*1e-9:e}")
    print("TargetVelocity corres. with Freq shift: (GHz) :",f"{TargetVelocityFreShift*1e-9:e}")
    return VPivot  

#tools to stop print on consol
def toggle_print(state=True):
    """Print active or not"""
    if state:
        sys.stdout = sys.__stdout__  # Activate  prints
    else:
        sys.stdout = open(os.devnull, 'w')  # desactivate Prints
    return    

"""
https://github.com/michael-betz/readTrc
Little helper class to load data from a .trc binary file.
This is the file format used by LeCroy oscilloscopes.
Thanks to M. Betz 09/2015
https://github.com/michael-betz/readTrc
"""

class Trc:
    _recTypes = (
        "single_sweep", "interleaved", "histogram", "graph",
        "filter_coefficient", "complex", "extrema",
        "sequence_obsolete", "centered_RIS", "peak_detect"
    )
    _processings = (
        "no_processing", "fir_filter", "interpolated", "sparsed",
        "autoscaled", "no_result", "rolling", "cumulative"
    )
    _timebases = (
        '1_ps/div', '2_ps/div', '5_ps/div', '10_ps/div', '20_ps/div',
        '50_ps/div', '100_ps/div', '200_ps/div', '500_ps/div', '1_ns/div',
        '2_ns/div', '5_ns/div', '10_ns/div', '20_ns/div', '50_ns/div',
        '100_ns/div', '200_ns/div', '500_ns/div', '1_us/div', '2_us/div',
        '5_us/div', '10_us/div', '20_us/div', '50_us/div', '100_us/div',
        '200_us/div', '500_us/div', '1_ms/div', '2_ms/div', '5_ms/div',
        '10_ms/div', '20_ms/div', '50_ms/div', '100_ms/div', '200_ms/div',
        '500_ms/div', '1_s/div', '2_s/div', '5_s/div', '10_s/div',
        '20_s/div', '50_s/div', '10Mouse left = +1 point  - Mouse right = -1 last point - Mouse middle = save&exit0_s/div', '200_s/div', '500_s/div',
        '1_ks/div', '2_ks/div', '5_ks/div', 'EXTERNAL'
    )
    _vCouplings = ('DC_50_Ohms', 'ground', 'DC_1MOhm', 'ground', 'AC,_1MOhm')
    _vGains = (
        '1_uV/div', '2_uV/div', '5_uV/div', '10_uV/div', '20_uV/div',
        '50_uV/div', '100_uV/div', '200_uV/div', '500_uV/div', '1_mV/div',
        '2_mV/div', '5_mV/div', '10_mV/div', '20_mV/div', '50_mV/div',
        '100_mV/div', '200_mV/div', '500_mV/div', '1_V/div', '2_V/div',
        '5_V/div', '10_V/div', '20_V/div', '50_V/div', '100_V/div',
        '200_V/div', '500_V/div', '1_kV/div'
    )

    def __init__(self):
        """
        use trc.open(fName) to open a Le Croy .trc file
        """
        self._f = None
        # offset to start of WAVEDESC block
        self._offs = 0
        self._smplFmt = "int16"
        self._endi = ""

    def open(self, fName):
        """
            _readS .trc binary files from LeCroy Oscilloscopes.
            Decoding is based on LECROY_2_3 template.
            [More info]
            (http://forums.ni.com/attachments/ni/60/4652/2/LeCroyWaveformTemplate_2_3.pdf)

            Parameters
            -----------
            fName = filename of the .trc file

            Returns
            -----------
            a tuple (x, y, d)

            x: array with sample times [s],

            y: array with sample  values [V],

            d: dictionary with metadata

            M. Betz 09/2015
        """
        print ("## Trc file extraction")
        with open(fName+".trc", "rb") as f:
            # Binary file handle
            self._f = f
            self._endi = ""
            temp = f.read(64)
            # offset to start of WAVEDESC block
            self._offs = temp.find(b'WAVEDESC')

            # -------------------------------
            #  Read WAVEDESC block
            # -------------------------------
            # Template name
            self._TEMPLATE_NAME = self._readS("16s", 16)
            if self._TEMPLATE_NAME != "LECROY_2_3":
                print(
                    "Warning, unsupported file template:",
                    self._TEMPLATE_NAME,
                    "... trying anyway"
                )
            # 16 or 8 bit sample format?
            if self._readX('H', 32):
                self._smplFmt = "int16"
            else:
                self._smplFmt = "int8"
            # Endian-ness ("<" or ">")
            if self._readX('H', 34):
                self._endi = "<"
            else:
                self._endi = ">"
            #  Get length of blocks and arrays
            self._lWAVE_DESCRIPTOR = self._readX("l", 36)
            self._lUSER_TEXT = self._readX("l", 40)
            self._lTRIGTIME_ARRAY = self._readX("l", 48)
            self._lRIS_TIME_ARRAY = self._readX("l", 52)
            self._lWAVE_ARRAY_1 = self._readX("l", 60)
            self._lWAVE_ARRAY_2 = self._readX("l", 64)

            d = dict()  # Will store all the extracted Metadata

            # ------------------------
            #  Get Instrument info
            # ------------------------
            d["INSTRUMENT_NAME"] = self._readS("16s", 76)
            d["INSTRUMENT_NUMBER"] = self._readX("l", 92)
            d["TRACE_LABEL"] = self._readS("16s", 96)

            # ------------------------
            #  Get Waveform info
            # ------------------------
            d["WAVE_ARRAY_COUNT"] = self._readX("l", 116)
            d["PNTS_PER_SCREEN"] = self._readX("l", 120)
            d["FIRST_VALID_PNT"] = self._readX("l", 124)
            d["LAST_VALID_PNT"] = self._readX("l", 128)
            d["FIRST_POINT"] = self._readX("l", 132)
            d["SPARSING_FACTOR"] = self._readX("l", 136)
            d["SEGMENT_INDEX"] = self._readX("l", 140)
            d["SUBARRAY_COUNT"] = self._readX("l", 144)
            d["SWEEPS_PER_ACQ"] = self._readX("l", 148)
            d["POINTS_PER_PAIR"] = self._readX("h", 152)
            d["PAIR_OFFSET"] = self._readX("h", 154)
            d["VERTICAL_GAIN"] = self._readX("f", 156)
            d["VERTICAL_OFFSET"] = self._readX("f", 160)
            # to get floating values from raw data:
            # VERTICAL_GAIN * data - VERTICAL_OFFSET
            d["MAX_VALUE"] = self._readX("f", 164)
            d["MIN_VALUE"] = self._readX("f", 168)
            d["NOMINAL_BITS"] = self._readX("h", 172)
            d["NOM_SUBARRAY_COUNT"] = self._readX("h", 174)
            # sampling interval for time domain waveforms
            d["HORIZ_INTERVAL"] = self._readX("f", 176)
            # trigger offset for the first sweep of the trigger,
            # seconds between the trigger and the first data point
            d["HORIZ_OFFSET"] = self._readX("d", 180)
            d["PIXEL_OFFSET"] = self._readX("d", 188)
            d["VERTUNIT"] = self._readS("48s", 196)
            d["HORUNIT"] = self._readS("48s", 244)
            d["HORIZ_UNCERTAINTY"] = self._readX("f", 292)
            d["TRIGGER_TIME"] = self._getTimeStamp(296)
            d["ACQ_DURATION"] = self._readX("f", 312)
            d["RECORD_TYPE"] = Trc._recTypes[
                self._readX("H", 316)
            ]
            d["PROCESSING_DONE"] = Trc._processings[
                self._readX("H", 318)
            ]
            d["RIS_SWEEPS"] = self._readX("h", 322)
            d["TIMEBASE"] = Trc._timebases[self._readX("H", 324)]
            d["VERT_COUPLING"] = Trc._vCouplings[
                self._readX("H", 326)
            ]
            d["PROBE_ATT"] = self._readX("f", 328)
            d["FIXED_VERT_GAIN"] = Trc._vGains[
                self._readX("H", 332)
            ]
            d["BANDWIDTH_LIMIT"] = bool(self._readX("H", 334))
            d["VERTICAL_VERNIER"] = self._readX("f", 336)
            d["ACQ_VERT_OFFSET"] = self._readX("f", 340)
            d["WAVE_SOURCE"] = self._readX("H", 344)
            d["USER_TEXT"] = self._readS(
                "{0}s".format(self._lUSER_TEXT),
                self._lWAVE_DESCRIPTOR
            )

            y = self._readSamples()
            y = d["VERTICAL_GAIN"] * y - d["VERTICAL_OFFSET"]
            x = np.arange(1, len(y) + 1, dtype=float)
            x *= d["HORIZ_INTERVAL"]
            x += d["HORIZ_OFFSET"]
        self.f = None
        self.Time = x
        self.Tension = y
        self.ScopeStatus = d
        #save data set in .csv"
        print("Save Tension(Time) .csv")
        np.savetxt(fName+".csv", np.vstack((x ,y)).T, delimiter=',')
        #np.savetxt(self.FName+"VelocityProfile.csv", np.vstack((self.VelocityProfile[:, 0] ,self.VelocityProfile[:, 1])).T, delimiter=',')
        
        return

    def _readX(self, fmt, adr=None):
        """ extract a byte / word / float / double from the binary file f """
        fmt = self._endi + fmt
        nBytes = struct.calcsize(fmt)
        if adr is not None:
            self._f.seek(adr + self._offs)
        s = struct.unpack(fmt, self._f.read(nBytes))
        if(type(s) == tuple):
            return s[0]
        else:
            return s

    def _readS(self, fmt="16s", adr=None):
        """ read (and decode) a fixed length string """
        temp = self._readX(fmt, adr).split(b'\x00')[0]
        return temp.decode()

    def _readSamples(self):
        # ------------------------
        #  Get main sample data with the help of numpys .fromfile(
        # ------------------------
        # Seek to WAVE_ARRAY_1
        self._f.seek(
            self._offs + self._lWAVE_DESCRIPTOR +
            self._lUSER_TEXT + self._lTRIGTIME_ARRAY +
            self._lRIS_TIME_ARRAY
        )
        y = np.fromfile(self._f, self._smplFmt, self._lWAVE_ARRAY_1)
        if self._endi == ">":
            y.byteswap(True)
        return y

    def _getTimeStamp(self, adr):
        """ extract a timestamp from the binary file """
        s = self._readX("d", adr)
        m = self._readX("b")
        h = self._readX("b")
        D = self._readX("b")
        M = self._readX("b")
        Y = self._readX("h")
        trigTs = datetime.datetime(
            Y, M, D, h, m, int(s), int((s - int(s)) * 1e6)
        )
        return trigTs

#!/usr/bin/env python
# coding: utf-8

# In[1]:
#get_ipython().run_line_magic('matplotlib', '')
#get_ipython().run_line_magic('matplotlib', 'inline')
#get_ipython().run_line_magic('matplotlib', 'tk')

from math import exp,pi,cos, sin
import numpy as np
import scipy
from scipy.signal import stft

from pylab import *

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import MultiCursor
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

from PyQt5.QtWidgets import QApplication, QMainWindow, QSlider, QLabel, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt

import csv
from pylab import *
import csv
import pandas as pd

import time

import os
import sys

import tkinter as tks
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

#matplotlib.use("TkAgg")
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
    def __init__(self,LambdaLaser,ChainResponse,Shift,FName,ShotNumber,nperseg):
        #all are calculate in SI but on application print in Nm and Ghz
        self.LambdaLaser=LambdaLaser*1e-9 #>SI
        self.Shift=Shift*1e9 #Hz PDV shift of reference Laser
        self.FName=FName #ShotName
        self.ChainResponse=ChainResponse*1e9
        self.ShotNumber=ShotNumber
        self.PDVFactor=self.LambdaLaser/2
        self.nperseg=nperseg
        self.STFTPDVWindow='hamming'
        
        
        #print ("Design*****")
        self.VPivot=self.LambdaLaser*self.Shift/2
        #print("VPivot (m/s) :",self.VPivot)
        self.MaxVelocityForChainResponse=self.PDVFactor*self.ChainResponse
        
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
                self.Time=np.append(self.Time,ti)
                self.Tension=np.append(self.Tension,vi)
                
    # Convertir les listes en arrays numpy
        # Convert lists to numpy arrays
        #self.Time = np.array(self.Time)
        #self.Tension = np.array(self.Tension)
    
        print("Number of points from DataLoad:", len(self.Time))
        return
    
    def PDVSetFrAcquisition(self):
        #Data extraction of acquisition sample rate in Sample/s
        print ("## PDVSetFrAcquisition calculation")
        self.Dtime=self.Time[2]-self.Time[1]
        self.FAcquisiton=1/self.Dtime
        print ("DTime (ns) :", self.Dtime*1e9)
        print ("FAquisition (GS/s) :",f"{self.FAcquisiton*1e-9:e}")
        return
    
    
    def SetPDVFFT(self):
        #Calculate FFT
        print ("## FTT signal calculation")
        self.N = len(self.Time)  # Taille du signal
        print ("Number of points: ", self.N)
        self.HSignalFFT= np.fft.rfft(self.Tension)  # Calcul de la FFT
        self.PDVSignalFFTTime= np.fft.rfftfreq(self.N, 1/self.FAcquisiton)  # Axe fréquentiel
        return
    
    def SetSTFTPDV(self,nperseg):
        # Calculate STFT (from segment nperseg (number of point), Sample rates,  output: Fefrequency : FPDV and related Time 
        #print ("## STFT signal calculation")
        self.SegSize=nperseg
        #print("neperseg: ",self.nperseg)
        #print("Window : "+self.STFTPDVWindow)
        self.WindowsSize=self.Dtime*self.SegSize
        #print (' Window (ns): ', self.WindowsSize*1e9) 
        self.FePDV, self.Time_stft, self.PDVSpectrogram = stft(self.Tension, self.FAcquisiton, nperseg=nperseg,window=self.STFTPDVWindow)
        return
    
    def SetVelocity(self):
        #Calculate PDVFactor
        print ("## Velocity signal calculation")
        self.Velocity=self.PDVFactor*self.FePDV
        return
   
    def NotebookGraph(self):   
        self.root = tk.Tk()
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(family="Verdana", size=12)
        self.root.title("PDV Analysis")
        self.shot_var = tk.StringVar()
        self.fname_var= tk.StringVar()
        self.nperseg_var=tk.StringVar()
        self.shot_dir=tk.StringVar()    
        self.ChainResponse_var=tk.StringVar()
        
        # Onglets
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=1)
    
        # Onglet Inputs at 
        self.frame_inputs = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_inputs, text="Data Load & Operation")
        self.CreateInputsTab(self.frame_inputs)
        
        #Console output
        self.CreateConsoleTab()
        sys.stdout = RedirectConsole(self.text_console)
        sys.stderr = RedirectConsole(self.text_console)
        
    
        # Ferme toute figure matplotlib résiduelle
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
    
    
    def CreateInputsTab(self, parent):
        style = ttk.Style()
        style.configure('TButton', font=('Arial', 12, 'bold'))
        # Variables
        self.shot_dir = tk.StringVar()
        self.fname = tk.StringVar()
        self.nperseg_var = tk.IntVar(value=500)
        self.ChainResponse_var = tk.DoubleVar(value=self.ChainResponse * 1e-9)
        Wavelength=self.LambdaLaser*1e9
        self.LambdaLaser_var= tk.StringVar(value=f'{Wavelength:.3f}')
        self.Shift_var= tk.DoubleVar(value=self.Shift * 1e-9)
        self.MaxVelocityForChainResponse_var=tk.StringVar(value=f'{self.MaxVelocityForChainResponse:.3f}')
        self.VPivot_var=tk.StringVar(value=f'{self.VPivot:.3f}')
        
        ttk.Label(parent, text="PDV Analysis", font=("Arial", 14, "bold")).pack(anchor="w", padx=10)
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=5,ipady=3)
        
        # Choix du répertoire (ShotNumber)
        ttk.Label(parent, text="Shot Directory:").pack(pady=5)
        frame_dir = tk.Frame(parent)
        frame_dir.pack()
        tk.Entry(frame_dir, textvariable=self.shot_dir, width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame_dir, text="Select Directory", command=self.select_directory).pack(side=tk.LEFT)
    
        # Choix du fichier (FName)
        ttk.Label(parent, text="File to analyse .csv with template (Time(s),Tension(V))").pack(pady=5)
        frame_file = tk.Frame(parent)
        frame_file.pack()
        tk.Entry(frame_file, textvariable=self.fname, width=50).pack(side=tk.LEFT, padx=5)
        tk.Button(frame_file, text="Select File ", command=self.select_file).pack(side=tk.LEFT)
    
        # nperseg
        ttk.Label(parent, text="Initial nperseg (FFT window size in number of point) :").pack(pady=5)
        tk.Entry(parent, textvariable=self.nperseg_var).pack()
        
        # Launch analysis
        ttk.Label(parent, text="Spectrogram, Raw datas and velocity figures are saved in png format ").pack(pady=15)
        ttk.Label(parent, text="Velocity data set in .csv file in ShotNumber Directory").pack(pady=5)
        ttk.Button(parent, text="Load Data Set for analysis", style='TButton',command=self.launch_analysis).pack(pady=5)
        
        ttk.Label(parent, text="PDV parameters", font=("Arial", 14, "bold")).pack(anchor="w", padx=10)
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=5)
        
        # LaserPDV Wavelength
        ttk.Label(parent, text="Laser PDV Wavelength (nm) :").pack(pady=15)
        tk.Entry(parent, textvariable=self.LambdaLaser_var, width=15).pack()
                
        # Frame contenant la ligne complète
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
        
        ttk.Button(parent, text="PDV Parameters Calculations", style='TButton',command=self.PDVParameters).pack(pady=15)
        
        ttk.Separator(parent, orient="horizontal").pack(fill="x", padx=10, pady=10, ipady=3)
        
        ttk.Button(parent, text="Exit", style='TButton', command=lambda: os._exit(0)).pack(pady=10)

    
    def select_directory(self):
        dirname = fd.askdirectory(title="Select Shot Directory")
        if dirname:
            self.shot_dir.set(dirname)
    
    def select_file(self):
        initialdir = self.shot_dir.get() if self.shot_dir.get() else "."
        filename = fd.askopenfilename(title="Select raw datas file .csv (Time(s), Tension(V))", initialdir=initialdir)
        if filename:
            self.fname.set(os.path.basename(filename))
            self.selected_file_fullpath = filename  # Chemin complet si besoin plus tard
            
        
    def write(self, string):
        self.output.insert(tk.END, string)
        self.output.see(tk.END)  # scroll automatique à la fin

    def flush(self):
        pass
    
    def CreateConsoleTab(self):
        self.frame_console = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_console, text="Console Output")
        self.text_console = tk.Text(self.frame_console, height=15, width=80)
        self.text_console.pack(fill='both', expand=True)
    
        sys.stdout = RedirectConsole(self.text_console)
        sys.stderr = RedirectConsole(self.text_console)
    
    def launch_analysis(self):
        
        for tab_id in self.notebook.tabs():
            tab_text = self.notebook.tab(tab_id, "text")
            if tab_text not in ("Data Load & Operation","Console Output"):
                self.notebook.forget(tab_id)
        
        
        # Récupère les valeurs saisies
        self.ShotNumber = self.shot_dir.get()
        self.FName = self.fname.get()
        self.nperseg = self.nperseg_var.get()
        
        #self.ChainResponse =self.ChainResponse_var.get()*1e9
        
        #print(self.ShotNumber,self.FName)
        # Charge tes données, calcule le STFT, etc.
        # (Tu insères ici SetSTFTPDV, LectureSignal, etc.)
        
        
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
        print("window type        :"+self.STFTPDVWindow)
        
        self.DataLoad(1)  # ← remplace par ta fonction réelle
        self.PDVSetFrAcquisition()
        self.SetPDVFFT()
        self.SetSTFTPDV(self.nperseg)
        
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
    
    
        # Onglet STFT Interactive
        self.frame_stft = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_stft, text="STFT Interactive")
        self.CreateSTFTPDVInteractive(self.frame_stft)
    
        # Onglet Raw Datas
        self.frame_graphs = ttk.Frame(self.notebook)
        self.notebook.add(self.frame_graphs, text="Raw Datas")
        self.NotebookGraphSpectrogram(self.frame_graphs)
    
        # Active l’onglet STFT directement
        self.notebook.select(self.frame_stft)

        for fig_num in plt.get_fignums():
            plt.close(fig_num)
            
    def CreateSTFTPDVInteractive(self, parent):
        self.param = self.nperseg
    
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2Tk(self.canvas, parent)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
    
        self.label = tk.Label(parent, text=f" STFT Window (nperseg) = {self.param} pt , {self.WindowsSize*1e9:.3f} ns")
        self.label.pack()
    
        self.ax.set_title("Spectrogram")
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Frequency (Hz)")
    
        self.STFTPDVWindow_var = tk.StringVar(value=self.STFTPDVWindow)
        fenetres = ['hann', 'hamming', 'blackman', 'bartlett', 'flattop']

        ttk.Label(parent, text="Windows STFT :").pack(anchor="w", padx=10, pady=(10, 0))
        combo = ttk.Combobox(parent, textvariable=self.STFTPDVWindow_var, values=fenetres, state="readonly")
        combo.pack(anchor="w", padx=10, pady=5)
        combo.bind('<<ComboboxSelected>>', lambda e: self.update_STFTPDVInteractiveplot(self.slider.get()))
    
        
        self.slider = ttk.Scale(parent, from_=64, to=2048, orient='horizontal')
        self.slider.set(self.param)
        self.slider.pack(fill=tk.X, padx=15, pady=15)
    
        self.slider.configure(command=self.update_STFTPDVInteractiveplot)
    
        self.Velocity_button = tk.Button(parent, text="ExtractVelocity", command=self.ExtractVelocityNotebook)
        self.Velocity_button.pack(pady=10)
    
        # Calcul initial
        self.SetSTFTPDV(self.param)
    
        # Affichage initial
        self.quadmesh = self.ax.pcolormesh(self.Time_stft, self.FePDV, np.abs(self.PDVSpectrogram), shading='gouraud')
        self.ax.set_ylim(min(self.FePDV), max(self.FePDV))
        self.ax.set_title("Spectrogram")
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        self.canvas.draw_idle()
        
                  
    def update_STFTPDVInteractiveplot(self, val):
        val = float(val)
        window = self.STFTPDVWindow_var.get()
        self.STFTPDVWindow = window
        self.param = int(val)
        self.nperseg=self.param
        self.label.config(
            text=f"STFT Window (nperseg) = {self.param} pt, {self.WindowsSize*1e9:.3f} ns, Window: {self.STFTPDVWindow}"
        )
    
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
    
        # Recalculer STFT
        self.SetSTFTPDV(self.param)
    
        # Effacer seulement le contenu des axes
        self.ax.clear()
    
        # Tracer à nouveau
        self.quadmesh = self.ax.pcolormesh(
            self.Time_stft,
            self.FePDV,
            np.abs(self.PDVSpectrogram),
            shading='gouraud'
        )
    
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Frequency (Hz)")
        self.ax.set_title("Spectrogram")
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.fig.savefig(self.FName+'Spectrogram.png')
        # Mettre à jour la figure dans Tkinter
        self.canvas.draw_idle()
    
    def ExtractVelocityNotebook(self):
        self.VelocityProfile = []
    
        # Fonction de clic dans le spectrogramme interactif
        def onclick(event):
            if event.inaxes == self.ax:  # Clic dans la figure STFT Interactive
                self.VelocityProfile.append((event.xdata, event.ydata*self.PDVFactor))
                #print(f"Point ajouté : {event.xdata:.4f}, {event.ydata:.4f}")
                self.ax.plot(event.xdata, event.ydata, 'rx')
                self.canvas.draw_idle()
    
        # Fonction de fin d'enregistrement : créer un nouvel onglet
        def stop_recording():
            self.canvas.mpl_disconnect(self.cid)
            print("Extraction is over")
            #print("Points extraits :", self.VelocityProfile)
            self.stop_button.destroy()
            # === Créer un nouvel onglet pour afficher les points ===
            self.frame_velocity = ttk.Frame(self.notebook)
            self.notebook.add(self.frame_velocity, text="Velocity Extraction")
            self.notebook.select(self.frame_velocity)
    
            # Figure vide pour affichage des points extraits
            fig_vel, ax_vel = plt.subplots(figsize=(6, 4))
            canvas_vel = FigureCanvasTkAgg(fig_vel, master=self.frame_velocity)
            canvas_vel.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            toolbar = NavigationToolbar2Tk(canvas_vel, self.frame_velocity)
            toolbar.update()
            toolbar.pack(side=tk.TOP, fill=tk.X)
            for fig_num in plt.get_fignums():
                plt.close(fig_num)
    
            ax_vel.set_title("Velocity profile")
            ax_vel.set_xlabel("Temps (s)")
            ax_vel.set_ylabel("Velocity (m/s)")
            ax_vel.grid(True)
    
            # Tracer les points extraits
            if self.VelocityProfile:
                x, y = zip(*self.VelocityProfile)
                ax_vel.plot(x, y, 'rx-', label="Velocityprofile Shot "+ self.FName)
                ax_vel.legend()
                print ('Save velocity profile in '+self.FName+"VelocityProfile.csv")
                np.savetxt(self.FName+"VelocityProfile.csv", np.vstack((x ,y)).T, delimiter=',')
                fig_vel.savefig(self.FName+'Velocity.png')
    
            canvas_vel.draw_idle()
    
        # Connexion du clic
        self.cid = self.canvas.mpl_connect('button_press_event', onclick)
    
        # Ajouter un bouton pour arrêter
        self.stop_button = tk.Button(self.root, text="End of Extraction", command=stop_recording)
        self.stop_button.pack(pady=10)
    
        #print("Cliquez sur le spectrogramme pour sélectionner des points.")
    
    def NotebookGraphSpectrogram(self, parent):
        fig, axs = plt.subplots(1, 2, sharex=False, sharey=False, figsize=(6, 4))
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)
        self.toolbar = NavigationToolbar2Tk(canvas, parent)
        self.toolbar.update()
        self.toolbar.pack(side=tk.TOP, fill=tk.X)
    
        axs[0].plot(self.Time, self.Tension)
        axs[0].set_title("Signal (t) " + self.FName)
        axs[0].set_xlabel("Time (s)")
        axs[0].set_ylabel("Amplitude (V)")
        axs[0].grid()
    
        axs[1].plot(self.PDVSignalFFTTime, np.abs(self.HSignalFFT))
        axs[1].set_title("FFT")
        axs[1].set_xlabel("Fe(Hz)")
        axs[1].set_ylabel("Magnitude")
        axs[1].set_yscale('log')
        axs[1].set_xlim(0, self.ChainResponse*1.5)
        axs[1].grid()
    
        fig.tight_layout()
        fig.savefig(self.FName+'RawData.png')
        for fig_num in plt.get_fignums():
            plt.close(fig_num)
        canvas.draw_idle()
    
    
    def runSTFTPDVInteractive(self):
        self.root.mainloop()    
      
    def PDVParameters(self):
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
        c.drawString(10, h - 61, "FAquisition (GS/s)   :"+str(self.FAcquisiton*1e-9)) 
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


def toggle_print(state=True):
    """Print active or not"""
    if state:
        sys.stdout = sys.__stdout__  # Activate les prints
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
    




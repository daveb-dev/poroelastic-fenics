#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 15:30:38 2019

@author: ranjeet
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_compare(file1,file2): 
    df1 = pd.read_csv(file1) 
    df2 = pd.read_csv(file2) 
    cols1 = df1.columns 
    cols2 = df2.columns 
    fig, ax = plt.subplots() 
    ax.plot(df1[cols1[2]],df1[cols1[0]], label="Analytical") 
    ax.plot(df2[cols2[2]],df2[cols2[0]],'o', label="Numerical")  
    ax.legend() 
    minor_ticks = np.arange(0.1,0.95,0.2) 
    minor_yticks = np.arange(0.1,1.2,0.2) 
    ax.set_xticks(minor_ticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True)  
    ax.grid(which='minor', alpha=0.5, linestyle='--') 
    plt.grid(True, linestyle='--') 
    plt.xlim(xmin=0) 
    plt.ylim(ymin=0, ymax=1.2) 
    plt.xlabel("Normalized horizontal distance, x/a") 
    plt.ylabel("Normalized pressure, p/p0") 
    plt.show() 


def plot_compare2(file1,file2, label1="A", label2="N"): 
    df1 = pd.read_csv(file1) 
    df2 = pd.read_csv(file2) 
    cols1 = df1.columns 
    cols2 = df2.columns  
    ax.plot(df1[cols1[2]],df1[cols1[0]], label=label1) 
    ax.plot(df2[cols2[2]],df2[cols2[0]],'o', label=label2)  

times= [10,700,1800, 3600,8400]
coords = [(0.85,1.05),(0.57,0.82), (0.35,0.8), (0.2,0.62),(0.2,0.25)] 


def plot(times, coords):
    fig, ax = plt.subplots()
    shapes =['o','s','p','^','*', 'D','d']
    colors = ['m','g','b','r','k']
    i=0
    
    #for file in Path('.').rglob('*a.csv'):
    for time in times:
        file1 = "test"+str(time)+"a.csv"  # Analytical
        file2 = "test"+str(time)+".csv"  # Numerical
        df1 = pd.read_csv(file1)
        #file2 = str(file)[:-5] +".csv"
        df2 = pd.read_csv(file2)
        cols1 = df1.columns 
        cols2 = df2.columns
        ax.plot(df1[cols1[2]],df1[cols1[0]], colors[i]) 
        ax.plot(df2[cols2[2]],df2[cols2[0]], colors[i]+'o')
        ax.text(coords[i][0], coords[i][1] ,"t = "+ str(time)+"s", fontsize=12, color=colors[i])
        i = i+1
    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    minor_ticks = np.arange(0.1,0.95,0.2) 
    minor_yticks = np.arange(0.1,1.2,0.2) 
    ax.set_xticks(minor_ticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True)  
    ax.grid(which='minor', alpha=0.5, linestyle='--') 
    plt.grid(True, linestyle='--') 
    plt.xlim(xmin=0) 
    plt.ylim(ymin=0, ymax=1.2) 
    plt.xlabel("Normalized horizontal distance, x/a",fontsize=12) 
    plt.ylabel("Normalized pressure, p/p0",fontsize=12) 
    
    
    
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
#    plt.xlabel(r"$\bar{x}$", fontsize=20)
#    plt.ylabel(r"$\bar{p}$", fontsize=20)
    
    #Create legend from custom artist/label lists
    
    #Create legend from custom artist/label lists
    #ax.legend([handle for i,handle in enumerate(handles) if i in display]+[simArtist,anyArtist],
    #          [label for i,label in enumerate(labels) if i in display]+['Simulation', 'Analytic'])
    
    plt.legend(lines,["Simulation", "Analytical"], fontsize=12)
    
    plt.show()
        
    
    

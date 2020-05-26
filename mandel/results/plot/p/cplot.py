import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def plotu(times, coords, keyword="p", xlabel="Pressure[Pa]", ylabel="y[m]"):
    fig, ax = plt.subplots()
    shapes =['o','s','p','^','*', 'D','d']
    colors = ['m','g','b','r','k']
    
    #for file in Path('.').rglob('*a.csv'):
    for i in range(len(times)):
        file1 = keyword+"_ana_"+str(i)+".csv"  # Analytical
        file2 = keyword+"_"+str(i)+".csv"  # Numerical
        df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        cols1 = df1.columns 
        cols2 = df2.columns
        ax.plot(df1.arc_length, df1[cols1[0]],  colors[i]) 
        ax.plot(df2.arc_length, df2[cols2[0]], colors[i]+'o')
        text = "t = "+times[i]
        ax.text(coords[i][0], coords[i][1] , text, fontsize=12, color=colors[i])
    

    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    minor_xticks = np.arange(0.1,0.95,0.2) 
    minor_yticks = np.arange(0,1.1, 0.1) 
    ax.set_xticks(minor_xticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True) 
    ax.grid(True)
    ax.grid(which='minor', alpha=0.5, linestyle='--') 
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    plt.xlim(xmin=0) 
    #plt.ylim(ymin=0, ymax=10.2) 
    plt.xlabel(xlabel,fontsize=14) 
    plt.ylabel(ylabel,fontsize=14) 
    
    plt.legend(lines,["Numerical solution", "Analytical solution"], fontsize=12)
    
    plt.show()
        

if __name__ == "__main__":
    xlabel = "Normalized pore pressure, p/p0"
    ylabel = "Normalized horizontal displacement, x/a"
    
    times=["1 min","10 min","50 min","100 min","200 min"]
    coords = [(0.8,1.0),(0.6,0.8),(0.4,0.6),(0.2,0.38),(0.2,0.15)] 
    plotu(times,coords, keyword="p", xlabel=xlabel, ylabel=ylabel)
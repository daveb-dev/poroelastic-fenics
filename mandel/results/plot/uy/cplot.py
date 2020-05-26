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
        ax.plot(df1.arc_length, df1[cols1[1]],  colors[i]) 
        ax.plot(df2.arc_length, df2[cols2[1]], colors[i]+'o')
        text = "t="+times[i]
        ax.text(coords[i][0], coords[i][1] , text, fontsize=12, color=colors[i])
    

    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    minor_xticks = np.arange(0.1,0.95,0.2) 
    minor_yticks = np.arange(-0.0001,0,0.00001) 
    ax.set_xticks(minor_xticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True) 
    ax.grid(True)
    ax.grid(which='both', alpha=0.5, linestyle='--') 
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    plt.xlim(xmin=0, xmax=1.25) 
    ax.ticklabel_format(axis="y",style='sci', scilimits=(0,0))
    #plt.ylim(ymin=0, ymax=10.2) 
    plt.xlabel(xlabel,fontsize=12) 
    plt.ylabel(ylabel,fontsize=12) 
    
    plt.legend(lines,["Numerical solution", "Analytical solution"], fontsize=12)
    
    plt.show()
        

if __name__ == "__main__":
    xlabel = "Normalized horizontal distance, x/b"
    ylabel = "Normalized vertical displacement, uy/b"
    
    times=["1min","10min","50min","100min","200min"]
    coords = [(0.9,-0.000056),(1.015,-0.000072),(1.02,-0.0000795),(1.02,-0.0000875),(0.9,-0.000098)] 
    plotu(times,coords, keyword="uy", xlabel=xlabel, ylabel=ylabel)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



#times= [10,700,1800, 3600,8400]
#coords = [(0.85,1.05),(0.57,0.82), (0.35,0.8), (0.2,0.62),(0.2,0.25)] 



def plot(times, coords, xlabel="Pressure[Pa]", ylabel="z[m]"):
    fig, ax = plt.subplots()
    shapes =['o','s','p','^','*', 'D','d']
    colors = ['m','g','b','r','k']
    i=0
    
    #for file in Path('.').rglob('*a.csv'):
    for time in times:
        #file1 = "test"+str(time)+"a.csv"  # Analytical
        file2 = "test"+str(time)+".csv"  # Numerical
        #df1 = pd.read_csv(file1)
        df2 = pd.read_csv(file2)
        #cols1 = df1.columns 
        cols2 = df2.columns
        #ax.plot(df1[cols1[0]], df1.arc_length, colors[i]) 
        ax.plot(df2[cols2[0]],df2.arc_length, colors[i]+'o')
        #ax.plot(df2[cols2[4]],df2[cols2[0]], colors[i]+'o')
        ax.text(coords[i][0], coords[i][0] ,"t = "+ str(time)+"s", fontsize=12, color=colors[i])
        i = i+1
    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    #minor_ticks = np.arange(0.1,0.95,0.2) 
    #minor_yticks = np.arange(0.1,1.2,0.2) 
    #ax.set_xticks(minor_ticks, minor=True) 
    #ax.set_yticks(minor_yticks, minor=True)  
    #ax.grid(which='minor', alpha=0.5, linestyle='--') 
    plt.grid(True, linestyle='--') 
    plt.xlim(xmin=0) 
    plt.ylim(ymin=0,ymax=10.0) 
    plt.xlabel(xlabel,fontsize=12) 
    plt.ylabel(ylabel,fontsize=12) 
    
    
    
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
        

if __name__ == "__main__":
    times= [0,1,2,3]
    coords = [(1.0,4.5e-5),(0,0),(0,0),(0,0)] 
    plot(times,coords)

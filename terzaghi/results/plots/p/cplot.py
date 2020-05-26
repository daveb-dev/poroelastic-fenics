import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



#times= [10,700,1800, 3600,8400]
#coords = [(0.85,1.05),(0.57,0.82), (0.35,0.8), (0.2,0.62),(0.2,0.25)] 



def plotp(times, coords, keyword="p", xlabel="Pressure[Pa]", ylabel="y[m]"):
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
        ax.plot(df1[cols1[0]], df1.arc_length, colors[i]) 
        ax.plot(df2[cols2[0]],df2.arc_length, colors[i]+'o')
        text = "t = "+times[i]
        print("text=", text)
        ax.text(coords[i][0], coords[i][1] , text, fontsize=12, color=colors[i])
    

    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    minor_xticks = np.arange(0,11000,1000) 
    minor_yticks = np.arange(0,10,1) 
    ax.set_xticks(minor_xticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True)  
    ax.grid(which='minor', alpha=0.5, linestyle='--') 
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    plt.xlim(xmin=0) 
    plt.ylim(ymin=0, ymax=10.2) 
    plt.xlabel(xlabel,fontsize=14) 
    plt.ylabel(ylabel,fontsize=14) 
    
    
    
#    ax.tick_params(axis='both', which='major', labelsize=12)
#    ax.tick_params(axis='both', which='minor', labelsize=12)
#    plt.xlabel(r"$\bar{x}$", fontsize=20)
#    plt.ylabel(r"$\bar{p}$", fontsize=20)
    
    #Create legend from custom artist/label lists
    
    #Create legend from custom artist/label lists
    #ax.legend([handle for i,handle in enumerate(handles) if i in display]+[simArtist,anyArtist],
    #          [label for i,label in enumerate(labels) if i in display]+['Simulation', 'Analytic'])
    
    plt.legend(lines,["Numerical solution", "Analytical solution"], fontsize=12)
    
    plt.show()
        

if __name__ == "__main__":
    
    times=["0.1 day","1 day", "3 day", "5 day", "10 day"]
    coords = [(6000,9.0),(6000,6.5),(5000,4),(3500,2),(1000,3)] 
    plotp(times,coords)

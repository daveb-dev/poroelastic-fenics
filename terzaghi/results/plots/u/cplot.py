import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



#times= [10,700,1800, 3600,8400]
#coords = [(0.85,1.05),(0.57,0.82), (0.35,0.8), (0.2,0.62),(0.2,0.25)] 



def plotu(times, coords, keyword="u", xlabel="Vertical displacement, uy [m]", ylabel="y[m]"):
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
        ax.plot(df1[cols1[1]], df1.arc_length, colors[i]) 
        ax.plot(df2[cols2[1]], df2.arc_length, colors[i]+'o')
        text = "t = "+times[i]
        ax.text(coords[i][0], coords[i][1] , text, fontsize=12, color=colors[i])
    

    
    #Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0,1,2)
    
    
    
    #Create custom artists
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='o', linestyle='')
    anyArtist = plt.Line2D((0,1),(0,0), color='k')
    lines=[simArtist,anyArtist]
    
    minor_xticks = np.arange(-0.0008,0,0.0001) 
    minor_yticks = np.arange(0,10,1) 
    ax.set_xticks(minor_xticks, minor=True) 
    ax.set_yticks(minor_yticks, minor=True)  
    ax.grid(which='minor', alpha=0.5, linestyle='--') 
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    plt.xlim(xmin=-.0008) 
    plt.ylim(ymin=0) 
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
    #coords = [(-0.006,9.0),(-0.002,6.5),(-0.002,4),(-0.006,2),(-0.004,3)] 
    coords = [(-0.0001,10),(-0.00025,9.3),(-0.00045,9.2),(-0.00062,9.8),(-0.0007,7.5)] 
    plotu(times,coords)

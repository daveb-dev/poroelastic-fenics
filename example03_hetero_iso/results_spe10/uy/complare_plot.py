import pandas as pd
import matplotlib.pyplot as plt

#0 864  
files =[1,5,10]

fig, ax = plt.subplots()

text_pos = [(173, -0.101e-2),(173,-0.641e-3), (173,-0.201e-3)]
#i=0
for file_ in files:
    df = pd.read_csv("uy_"+str(file_)+".csv")
    cols = df.columns
    label="at " + str(file_*864)+"s"
    ax.plot(df[cols[4]], df[cols[1]], label=label)
    #plt.text(173, text_pos[i][1], s=label, color='red')
    #i = i + 1

    
#ax.set_xticks(minor_ticks, minor=True) 
#ax.set_yticks(minor_yticks, minor=True)  
#ax.grid(which='minor', alpha=0.5, linestyle='--') 
plt.grid(True, linestyle='--')
#plt.xlim(xmin=0)
#plt.ylim(ymin=0, ymax=1.2)
#plt.xlabel("Horizontal distance (m).", fontsize=12)
plt.xlabel("Horizontal distance (m)")
plt.ylabel("Vertical displacement (m)")
#plt.legend()

ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='y')
plt.text(173, -0.101e-2, 't=8640s')
plt.text(173, -0.641e-3, 't=4320s')
plt.text(173, -0.207e-3, 't=864s')

plt.show()

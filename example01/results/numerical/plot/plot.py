import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


fig, ax = plt.subplots()
colors= {'A':'r', 'B':'b', 'C':'g'}

files = ["uy_A", 'uy_B', 'uy_C']
linestyles = ["-","--", "-."]
i=0
for filepath in files:
    layer = filepath[-1]
    df1 = pd.read_csv(filepath+"_600s.csv")
    df2 = pd.read_csv(filepath+"_1800s.csv")
    label1 = "$u_y$ at " + layer + " after " + "10min"
    label2 = "$u_y$ at " + layer + " after " + "30min"
    ax.plot(df1.arc_length, df1['u:1'], ":"+colors[layer], label=label1)
    ax.plot(df2.arc_length, df2['u:1'], "-"+colors[layer], label=label2)
    i = i + 1



#plt.annotate('Test', xy=(-0.01,0.2), xytext=(-0.07, 0.3),
#            arrowprops=dict(facecolor='black', arrowstyle='-|>'))

#plt.text(1.0, -0.085, "A")
#plt.text(1.0, -0.0385, "A")
#plt.text(1.0, -0.048, "B")
#plt.text(1.0, -0.0252, "B")
#plt.text(1.0, -0.01, "C")
#plt.text(1.0, -0.0193, "C")
#
minor_xticks = np.arange(0.1, 0.95, 0.1)
minor_yticks = np.arange(-0.03, -0.14, -0.02) 
ax.set_xticks(minor_xticks, minor=True)
ax.set_yticks(minor_yticks, minor=True)
ax.grid(which='both', alpha=0.5, linestyle='--') 
#plt.gca().invert_yaxis()
plt.grid(True)
plt.xlabel("Distance along interface (m)", fontsize=12)
plt.ylabel("Vertical displacement (m)", fontsize=12)
plt.legend(["at 10min","at 30min"], ncol=2)


## Shrink current axis's height by 10% on the bottom
#box = ax.get_position()
#ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                 box.width, box.height * 0.9])
#
#
## Put a legend below current axis
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
#          fancybox=True,ncol=3,handleheight=1)
#

plt.annotate('A', xy=(0.9,-0.085), xytext=(1.0,-0.065),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))
plt.annotate('A', xy=(0.9,-0.038), xytext=(1.0,-0.065),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))

plt.annotate('B', xy=(0.2,-0.0538), xytext=(0.1,-0.045),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))
plt.annotate('B', xy=(0.2,-0.034), xytext=(0.1,-0.045),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))

plt.annotate('C', xy=(0.6,-0.0196), xytext=(0.52,-0.0168),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))
plt.annotate('C', xy=(0.6,-0.0106), xytext=(0.52,-0.0168),
            arrowprops=dict(facecolor='black', arrowstyle='-|>'))


#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.axis([0.0, None, -0.1, None])
plt.show()
#uy_1h_C.csv
#uy_2h_A.csv
#uy_1h_A.csv
#uy_2h_B.csv
#uy_1h_B.csv
#uy_2h_C.csv

if __name__ == "__main__":
    import os
    print(os.listdir('.'))
    
    

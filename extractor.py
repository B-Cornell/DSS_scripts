import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
outgoingvel = []
incomingvel = []
outgoingvelprob = []
incomingvelprob = []

#this unpacks the data file
#173568 pairs
pairid, massa, massb, separation, onvel, offvel, incoming, relvelangle, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31, a32, a33, a34, a35, a36, a37, a38, a39, a40, a41, a42, a43, a44, a45, a46, a47, a48, a49, a50, prob = np.loadtxt('pyfile.csv', delimiter = ',', unpack = True)

def vel(i):
    total = math.sqrt(onvel[i]**2 + offvel[i]**2)
    return(total)

for i in range(len(incoming)):
    if incoming[i] == 1:
 
        incomingvel.append(vel(i))
        incomingvelprob.append(vel(i)*prob[i])
        
    elif incoming[i] == 0:
      
        outgoingvel.append(vel(i))
        outgoingvelprob.append(vel(i)*prob[i])
        

plt.hist(outgoingvel,color = 'b',bins = 30,stacked = False,histtype = 'bar',rwidth = .7)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Counts')
plt.title('Outgoing')
plt.show()
plt.hist(outgoingvelprob,color = 'b',bins = 30,stacked = False,histtype = 'bar',rwidth = .7)
plt.xlabel('Velocity (km/s)')
plt.ylabel('Counts Weighted with Probability')
plt.title('Weighted Counts\nOutgoing')
plt.show()





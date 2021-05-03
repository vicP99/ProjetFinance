import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

Y=[0.0990432,0.128834,0.153075,0.155852]
Y_anti=[0.117645,0.138412,0.158107,0.153669]
X=[10,100,1000,10000]


X_u=[-0.00475319,0.094001,0.138155,0.151381]
X_u_anti=[0.0717875,0.119784,0.15085,0.151436]
X_o_anti=[0.163503,0.157039,0.165364,0.155902]
X_o=[0.20284,0.163668,0.167995,0.160323]

plt.plot(X,Y,label="estimation par monte carlo")
plt.xscale("log")
plt.xlabel("nombre de simulation")
plt.ylabel("P")
plt.axhline(y=0.156055,linestyle="--",color="green",label="P analytique")
for i in range(4):
    if(i==0):
        plt.plot([10**(i+1),10**(i+1)],[X_u[i],X_o[i]],linestyle="dotted",color="red",label="intervalle de confiance à 95%")
    else:
        plt.plot([10**(i+1),10**(i+1)],[X_u[i],X_o[i]],linestyle="dotted",color="red")
plt.plot(X,Y_anti,label="estimation antithétique")
for i in range(4):
    if(i==0):
        plt.plot([10**(i+1),10**(i+1)],[X_u_anti[i],X_o_anti[i]],linestyle="-",color="purple",label="intervalle de confiance à 95%")
    else:
        plt.plot([10**(i+1),10**(i+1)],[X_u_anti[i],X_o_anti[i]],linestyle="-",color="purple")
plt.legend()
plt.show()

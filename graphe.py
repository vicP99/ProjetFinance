import matplotlib.pyplot as plt
import matplotlib as cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.ticker import LinearLocator, FormatStrFormatter
""" 
fichier=open("out/sortie.txt","r")
Y=[]
E=[]
i=0
compteur=0
for ligne in fichier:
    if ligne!='' :
        if compteur==0:
            Y.append(float(ligne[:-1]))
            compteur=1
        else:
            E.append(float(ligne[:-1]))
            compteur=0
        i+=1

E=np.asarray(E)
Y=np.asarray(Y)
sup=Y+E
inf=Y-E
print(len(E))
X=np.linspace(-1,1,101)
plt.plot(X,Y,":",color="r",label="valeur")
plt.fill_between(X,inf,sup,color='#1E90FF',label="IC",alpha=0.75)
plt.legend()
plt.title("Tracé du payoff de l'option spread pour K dans [-1,1]")
plt.xlabel("rho")
plt.ylabel("Valeur de l'option (P)") 
plt.show()


 




 fichier=open("out/sortie.txt","r")
Val=[]
ValCond=[]
E=[]
Econd=[]
i=0
compteur=0
for ligne in fichier:
    if ligne!='' :
        if compteur==0:
            Val.append(float(ligne[:-1]))
            compteur+=1
        elif compteur==1:
            ValCond.append(float(ligne[:-1]))
            compteur+=1
        elif compteur==2:
            E.append(float(ligne[:-1]))
            compteur+=1
        elif compteur==3:
            Econd.append(float(ligne[:-1]))
            compteur=0

Val=np.asarray(Val)
ValCond=np.asarray(ValCond)
E=np.asarray(E)
Econd=np.asarray(Econd)

sup=Val+E
inf=Val-E

supCond=ValCond+Econd
infCond=ValCond-Econd"""
fig = plt.figure()
ax = fig.gca(projection='3d')
fichier=open("out/sortie.txt","r")
Val=[]
A=np.eye(101,101)
i=0
j=0
compteur=0
for ligne in fichier:
    if ligne!='' :
        if(j==101):
            j=0
            i+=1
        A[i,100-j]=float(ligne[:-1])
        j+=1


sigma1=np.linspace(0,0.8,101)
sigma2=np.linspace(0,0.8,101)
X,Y=np.meshgrid(sigma1,sigma2)
a=sigma1
b=A()
print(Y)
#plt.imshow(A)
surf = ax.plot_surface(X, Y, A)
#plt.colorbar()
#plt.xticks([0,1:2.5,25,37.5,50,62.5,75,87.5,100],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
#plt.yticks([100,87.5,75,62.5,50,37.5,25,12.5,0],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
plt.title("Prix de l'option spread en fonction de la variance")
plt.xlabel("sigma2")
plt.ylabel("sigma1") 
plt.show()  
"""

X=np.linspace(1000,100000,100)

P=0.156055
#plt.plot(X,E)
#plt.plot(X,Econd)
#plt.plot(X,E/Econd)
#plt.title("Erreur classique/ erreur avec conditionnement")
#plt.xlabel("Nombre de simulation")
plt.plot(X,P*np.ones(100),":",color="r",label="valeur réelle")
#plt.fill_between(X,inf,sup,color='#1E90FF',label="Monte Carlo classique")
plt.fill_between(X,infCond,supCond,color='#90EE90',label="Monte Carlo avec condtionnement",alpha=0.75)
plt.legend()
plt.title("Visualisation des intervalles de confiance à 90%")
plt.xlabel("Nombre de simulation")
plt.ylabel("Valeur de l'option (P)") 
plt.show() """
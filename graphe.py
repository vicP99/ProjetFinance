import matplotlib.pyplot as plt
import numpy as np

""" fichier=open("out/sortie.txt","r")
Y=[]
i=0
for ligne in fichier:
    if ligne!='' :
        Y.append(float(ligne[:-2]))
        i+=1

X=[]
j=0
for j in range (i):
    X.append(-1 + (j+1)*2/i)

plt.plot(X,Y)
plt.show()
 """
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

for i in range(100):
    if (E[i]<0):
        print(i)

Val=np.asarray(Val)
ValCond=np.asarray(ValCond)
E=np.asarray(E)
Econd=np.asarray(Econd)

sup=Val+E
inf=Val-E

supCond=ValCond+Econd
infCond=ValCond-Econd

X=np.linspace(1000,100000,100)

P=0.156055
#plt.plot(X,E)
#plt.plot(X,Econd)
plt.plot(X,E/Econd)
plt.title("Erreur classique/ erreur avec conditionnement")
plt.xlabel("Nombre de simulation")
""" plt.plot(X,P*np.ones(100),":",color="r",label="valeur réelle")
plt.fill_between(X,inf,sup,color='#1E90FF',label="Monte Carlo classique")
plt.fill_between(X,infCond,supCond,color='#90EE90',label="Monte Carlo avec condtionnement",alpha=0.75)
plt.legend()
plt.title("Visualisation des intervalles de confiance à 90%")
plt.xlabel("Nombre de simulation")
plt.ylabel("Valeur de l'option (P)") """
plt.show()
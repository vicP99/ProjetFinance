import matplotlib.pyplot as plt

fichier=open("out/sortie.txt","r")
Y=[]
i=0
for ligne in fichier:
    ligne=fichier.readline()
    if ligne!='' :
        Y.append(float(ligne[:-2]))
        i+=1

X=[]
j=0
for j in range (i):
    X.append(-1 + (j+1)*2/i)

plt.plot(X,Y)
plt.show()

close(fichier)



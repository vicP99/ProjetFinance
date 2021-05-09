import numpy as np
import math as m
import matplotlib.pyplot as plt

np.random.seed(1)
dt=0.01
sigma1=0.25
sigma2=0.3
r=0.01
rho=0.5

t=np.linspace(0,10,1001)
x=np.linspace(-10,10,1001)
N1=np.random.normal(0,m.sqrt(dt),1000)
N3=np.random.normal(0,m.sqrt(dt),1000)

W1=np.concatenate((np.array([0]),np.cumsum(N1)))
W3=np.concatenate((np.array([0]),np.cumsum(N3)))
W2=rho*W1 + m.sqrt(1-rho*rho)*W3

S1=np.exp((r-sigma1**2/2)*t + sigma1*W1)
S2=np.exp((r-sigma2**2/2)*t + sigma1*W2)

K=0.2

plt.plot(t,S1,color='r',label='S1')
plt.plot(t,S2,color='b',label='S2')
#np.maximum(S1,0)
plt.plot(t,np.maximum(S1,S2)-K,color='g',label='max(S1,S2)-K')
""" plt.plot(x,np.maximum(x-K,0),label="(x-K)+")
plt.plot(x,np.maximum(K-x,0),label="(K-x)+") """
plt.xlabel("t")
plt.title("Option BestOf")
plt.legend()
plt.show()
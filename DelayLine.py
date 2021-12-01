# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 09:57:46 2021

@author: Johan
"""

import numpy as np
import matplotlib.pyplot as plt


#%%######################COAXIAL CABLE #################################

#R L C linéique 
lR = 1e-3 #ohm/m
lC = 50e-12 #F/m
lL = 300e-9 #H/m
#http://rfelektronik.se/manuals/Datasheets/Coaxial_Cable_Attenuation_Chart.pdf

Lline = 100 #for specific line length calculation
Lline2 = np.linspace(1,1000,1000) #for plotting graph
w = np.linspace(1,1e9, 10000)

ltau = (lL*lC)**0.5 #s/m

w0 = 1/((Lline*lL*Lline*lC)**0.5) #résonnance pulsation
Q = (1/(lR*Lline))*((lL/lC)**0.5) #facteur de qualité

f0fLline = (1/((Lline2*lL*Lline2*lC)**0.5))/(2*np.pi) #cut off frequency as a function of line length
taufLline = (Lline2*lL*Lline2*lC)**0.5 #delay as a function of line length

H = 1/(1+1j*(w/(Q*w0))-((w/w0)**2)) #transfer function of delay line

gain = 10*np.log(np.abs(H)) #gain in dB
cut_off_freq = 0
loss_100MHz = 0
for k in range(len(w)):
    if gain[k]<-3 and cut_off_freq ==0:
        cut_off_freq = w[k]/(2*np.pi)
    if k>0 and k<len(w)-1:
        if w[k+1]/(2*np.pi)>100e6 and w[k-1]/(2*np.pi)<100e6 or w[k]/(2*np.pi)==100e6:
            loss_100MHz = gain[k]

plt.figure() #gain bode diagram of line
plt.semilogx(w/(2*np.pi),gain,label = 'gain without linear loss')
plt.xlabel('f (Hz)')
plt.ylabel('gain (dB)')
plt.title('Résonnance = '+str(np.around(w0/(2*np.pi),0))+' Hz')
plt.legend()

print('R total = '+str(lR*Lline)+' ohm')
print('L total = '+str(lL*Lline*10**3)+' mH')
print('C total = '+str(lC*Lline*10**9)+ ' nF')
print('Losses @100MHz = '+str(round(loss_100MHz,0)))
print('Cut off freq @-3dB= '+str(np.around(cut_off_freq,0)))
print('tau = '+str(ltau*Lline*10**9)+' ns ---> '+str(ltau*Lline*156e6)+' noeuds')

fig, ax1 = plt.subplots() # cut off frequency and delay as a function of line length
ax1.set_xlabel('L line')
ax1.set_ylabel('f0 (MHz)',color = 'red')
ax1.plot(Lline2,f0fLline,color='red')
ax2 = ax1.twinx()
ax2.set_xlabel('L line')
ax2.set_ylabel('tau',color = 'blue')
ax2.plot(Lline2,taufLline,color='blue')
fig.tight_layout()


#%%#################### Discrete line

R = 50
L = 2e-7
C = 3e-12

N_cell = 5000

tau = N_cell*(L*C)**0.5 #s

w0 = 1/((L*C)**0.5)
Q = (1/(R))*((L/C)**0.5)

H = 1/(1+1j*(w/(Q*w0))-((w/w0)**2))

gain = 20*np.log(np.abs(H))

plt.figure()
plt.semilogx(w/(2*np.pi),gain,label = 'gain')
plt.xlabel('f (Hz)')
plt.ylabel('gain (dB)')
plt.title('Résonnance = '+str((np.around(w0/(2*np.pi),0))*10**-6)+' MHz')
plt.legend()

print('R = '+str(R)+' ohm')
print('L = '+str(L*10**3)+' mH')
print('C = '+str(C*10**9)+ ' nF')
print('Cut off freq = '+str(np.around(w0/(2*np.pi),0)))
print('delay = '+str(tau*10**9)+' ns ---> '+str(np.around((tau*(156e6)),1))+' noeuds')

C2 = np.linspace(0.5e-12,1e-9,1000)
L2 = np.linspace(0.022e-6,10e-6,1000)
f0fC2 = ((1/((L*C2)**0.5))/(2*np.pi))*(10**-6) #MHz
f0fL2 = ((1/((L2*C)**0.5))/(2*np.pi))*(10**-6) #MHz
taufC2 = (N_cell*(L*C2)**0.5)*(10**9)*(10**-9)*(156e6) #ns *(10**-9)*(156e6) #for number of nodes
taufL2 = (N_cell*(L2*C)**0.5)*(10**9)*(10**-9)*(156e6) #ns *(10**-9)*(156e6) #for number of nodes

fig, ax1 = plt.subplots()
ax1.set_xlabel('C')
ax1.set_ylabel('f0 (MHz)',color = 'red')
ax1.plot(C2,f0fC2,color='red')
ax2 = ax1.twinx()
ax2.set_xlabel('C')
ax2.set_ylabel('noeuds',color = 'blue')
ax2.plot(C2,taufC2,color='blue')
fig.tight_layout()
plt.title('L = '+str(L)+', N_cell = '+str(N_cell))

fig, ax1 = plt.subplots()
ax1.set_xlabel('L')
ax1.set_ylabel('f0 (MHz)',color = 'red')
ax1.plot(L2,f0fL2,color='red')
ax2 = ax1.twinx()
ax2.set_xlabel('L')
ax2.set_ylabel('noeuds',color = 'blue')
ax2.plot(L2,taufL2,color='blue')
fig.tight_layout()
plt.title('C = '+str(C)+', N_cell = '+str(N_cell))
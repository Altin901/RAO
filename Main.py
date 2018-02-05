
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 15:57:02 2016

@author: noudd
"""


import numpy as np
import matplotlib.pyplot as plt
import math

import pandas as pd
import os
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.special as special
from scipy.stats import norm
from scipy.interpolate import interp1d



# Units (SI): m,s,rad,
Spectrum = "JONSWAP"	# JONSWAP or PM


'''
Hs1 = float(input('What is the significant wave height? '))
print(hs1)
'''


Hs = 1.75
Tp = 8*1.2859
angle = 22.5
γ = 3.3
xT = 11

COM = [21.75,0,7.64]
P1 = [21.75,11.25,40.20]

#   Amount of RAO to be imported
step = math.floor(360 / angle)  

#Import all files from python file location
path = os.getcwd()
path = os.chdir("DATA_INPUT")

files = os.listdir(path)
files

#Import only excel files from python file location
files_xls = [f for f in files if f[-3:] == 'xls']
files_xls

#Determine size of excel files [column, rows]
RAOn = pd.read_excel(files_xls[0])
RAO = np.zeros((RAOn.shape[0],RAOn.shape[1],step))
μ = np.zeros((step))

#Store excel files in 3d array
for x in range(0,step):
    RAO[:,:,x] = pd.read_excel(files_xls[x])
    μ[x] = x * angle


ω = RAO[:,0,0]
ω_p = 2*math.pi/Tp
A_γ = 1-0.287*math.log(γ)

σ = np.zeros(np.size(ω))
for n in range(0,np.size(ω)):
    σ[n] = 0.07 if ω[n] <= ω_p else 0.09




S_PM = 5/16 * Hs**2 * ω_p**4*ω**-5 * np.exp(-5/4 * (ω/ω_p)**-4)
S_j = A_γ * S_PM * γ**(np.exp(-0.5* ((ω - ω_p)/ (σ * ω_p) )**2))

if Spectrum == "JONSWAP":	S_ζ = S_j
if Spectrum == "PM" :		S_ζ = S_PM	

f0 = interp1d(ω, S_ζ, kind='cubic')
f1 = interp1d(ω, S_ζ*ω, kind='cubic')
f2 = interp1d(ω, S_ζ*ω**2, kind='cubic')







'''
xnew = np.linspace(0, 10, num=41, endpoint=True)
plt.plot(ω, S_ζ, 'o', ω, f0(ω), '--')
plt.legend(['data', 'cubic'], loc='best')
'''





m0 = integrate.quad(f0,ω[0],(ω[(len(ω))-1]))
m1 = integrate.quad(f1,ω[0],(ω[(len(ω))-1]))
m2 = integrate.quad(f2,ω[0],(ω[(len(ω))-1]))

ν = math.sqrt(((m0[0]*m2[0])/m1[0]**2)-1)
Tz = 2*math.pi*math.sqrt(m0[0]/m2[0])
T1 = 2*math.pi*m0[0]/m1[0]
C1 = 1+ (ν**2/(1+ν**2)**(3/2))
C2 = 1/2*ν/(1+ν**2)
μT = C1*T1


Hw = 2.15*Hs/2
σT = C2*Hs/Hw*T1
ρTH = norm.cdf((xT-μT)/σT)


Sr = np.zeros((np.size(ω),6,step))



for n in range(0,len(ω)):
	for m in range(0,step-1):
		Sr[n,:,m] = RAO[n,1,m]**2*S_ζ[n],RAO[n,3,m]**2*S_ζ[n],RAO[n,5,m]**2*S_ζ[n],RAO[n,7,m]**2*S_ζ[n],RAO[n,9,m]**2*S_ζ[n],RAO[n,11,m]**2*S_ζ[n]

print(Sr)

Srx = Sry = Srz = SrΦ = Srθ = SrΨ = np.zeros((np.size(ω),step))
m0x = m0y = m0z = m0Φ = m0θ = m0Ψ = np.zeros(step)




for m in range(0,step-1):
	for n in range(0,len(ω)):
		Srx[n,m] = RAO[n,1,m]**2*S_ζ[n]
		Sry[n,m] = RAO[n,3,m]**2*S_ζ[n]
		Srz[n,m] = RAO[n,5,m]**2*S_ζ[n]
		SrΦ[n,m] = RAO[n,7,m]**2*S_ζ[n]
		Srθ[n,m] = RAO[n,9,m]**2*S_ζ[n]
		SrΨ[n,m] = RAO[n,11,m]**2*S_ζ[n]
	f0x = interp1d(ω,Srx[:,m], kind='linear')	
	f0y = interp1d(ω,Sry[:,m], kind='linear')	
	f0z = interp1d(ω,Srz[:,m], kind='linear')	
	f0Φ = interp1d(ω,SrΦ[:,m], kind='linear')
	f0θ = interp1d(ω,Srθ[:,m], kind='linear')
	f0Ψ = interp1d(ω,SrΨ[:,m], kind='linear')	
	Inx = integrate.quad(f0x,ω[0],(ω[(len(ω))-1]))
	Iny = integrate.quad(f0y,ω[0],(ω[(len(ω))-1]))
	Inz = integrate.quad(f0z,ω[0],(ω[(len(ω))-1]))
	InΦ = integrate.quad(f0Φ,ω[0],(ω[(len(ω))-1]))
	Inθ = integrate.quad(f0θ,ω[0],(ω[(len(ω))-1]))
	InΨ = integrate.quad(f0Ψ,ω[0],(ω[(len(ω))-1]))
	m0x[m] = Inx[0]
	m0y[m] = Iny[0]
	m0z[m] = Inz[0]
	m0Φ[m] = InΦ[0]
	m0θ[m] = Inθ[0]
	m0Ψ[m] = InΨ[0]



x = 2.14 * 2* math.sqrt(np.amax(m0x))
y = 2.14 * 2* math.sqrt(np.amax(m0y))
z = 2.14 * 2* math.sqrt(np.amax(m0z))
Φ = 2.14 * 2* math.sqrt(np.amax(m0Φ))
θ = 2.14 * 2* math.sqrt(np.amax(m0θ))
Ψ = 2.14 * 2* math.sqrt(np.amax(m0Ψ))

print(x,y,z,Φ,θ,Ψ)


'''


plt.figure(1)
plt.plot((ω),S_PM,label="Bretschneider spectrum")
plt.plot((ω),S_j,label="Jonswap spectrum")
plt.legend(loc=0)
plt.savefig('wave spectrum.png')


plt.show()

'''
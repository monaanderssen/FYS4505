
from math import *
import numpy as np 
from matplotlib import pyplot as plt 

#Parameters for Silicon:
rho_Si = 2.33 #g/cm^3
Z_Si = 14.
A_Si = 28 # u = g/mol
I_Si = 173e-6#10e-6*Z_Si #MeV

#Parameters for Germanium
rho_Ge = 5.323 #g/cm^3
Z_Ge = 32.
A_Ge = 72.64 # u = g/mol
I_Ge = 350e-6#10e-6*Z_Ge #MeV

#Div parameters:
N_a = 6.022e23 #mol^-1
r_e = 2.82e-13 #cm
m_e = 0.511 #MeV/c^2
K = 2*pi*N_a*r_e**2*m_e #cm^2 MeV/mol

# Mass alpha particle:
M = 3.727e3 #MeV/c^2

#Mass proton:
m_p = 938.3 #MeV/c^2

def bethe(rho, Z, A, I, E, type):
	if type == 0:
		beta = sqrt(E*2/M)
		gamma = 1/sqrt(1-beta**2)
		T_max = 2.*m_e*(beta*gamma)**2/(1+(2*m_e)/M *gamma + m_e**2/M**2) #MeV
		z = 2.
	elif type ==1:
		gamma = E/m_e**2 +1
		beta = np.sqrt(1-1/gamma**2)
		T_max = 5.
		z = 1.
	elif type ==2:
		gamma = E/m_e + 1
		beta = sqrt(1 - gamma**(-2))
		T_max = 2.*m_e*(beta*gamma)**2/(1+(2*m_e*gamma)/m_p  + m_e**2/m_p**2) #MeV
		z= 1.
	S = K*Z*z**2/(A*(beta**2))*(np.log(2*m_e*(beta*gamma)**2*T_max/I**2) - 2*beta**2) #Kan ogsa gange inn rho, men da blir enhetene annerledes!
	return S

S_Si = bethe(rho_Si, Z_Si, A_Si, I_Si, E=5., type=0)
S_Ge = bethe(rho_Ge, Z_Ge, A_Ge, I_Ge, E=5., type=0)

S_Si_e = bethe(rho_Si, Z_Si, A_Si, I_Si, E=5., type=1)
S_Ge_e = bethe(rho_Ge, Z_Ge, A_Ge, I_Ge, E=5., type=1)

print 'The stopping power for an alpha particle in silicon is', S_Si, 'MeVcm^2/g'
print 'The stopping power for an alpha particle in germanium is', S_Ge, 'MeVcm^2/g'
print 'The stopping power for an electron in silicon is',S_Si_e, 'MeVcm^2/g'
print 'The stopping power for an electron in germanium is',S_Ge_e, 'MeVcm^2/g'

# Exercise 2

#Parameters for aluminum:
Z_Al = 13.
rho_Al = 2.6989 #g/cm^3
I_Al = 166e-6
A_Al = 26.98

#Parameters for Portland concrete:
Z_con = 11.099
rho_con = 2.3 #g/cm^3
I_con = 132.047e-6
A_con = 95.33

#Div parameters:
N = 100000
c = 2.99e8 #m/s

S_Al = np.zeros(N)
S_con = np.zeros(N)

for i in range(N):
	E = np.linspace(0.1,10e2,N)
	S_Al[i] = bethe(rho_Al, Z_Al, A_Al, I_Si, E[i], type=2)
	S_con[i] = bethe(rho_con, Z_con, A_con, I_con, E[i], type=2)

plt.semilogx(E, S_Al)
plt.semilogx(E, S_con)
plt.xlabel('Energy [MeV]')
plt.ylabel('Stopping power [Mev cm^2/g]')
plt.title('Energy loss of the proton')
plt.legend(['Aluminum', 'Portland concrete'])
plt.show()

## Range

from scipy.integrate import simps

def InverseStoppingPower(Energy, rho, I, A, Z):
	electronmass = 0.511 #MeV
	particlemass = 938.3 #MeV
	z = 1.
	K = 1e4*2*pi*6.022e23*2.82e-15**2*electronmass
	gamma = Energy/particlemass + 1
	electronmass = 0.511 #Mev
	betaSquare = 1-1/gamma**2
	Tmax = 2*electronmass*betaSquare*gamma**2/( 1 + 2*(electronmass/particlemass) + (electronmass/particlemass)**2 )
	S = rho*K*(Z/A)*z**2/betaSquare*( log(2*electronmass*betaSquare*gamma**2*Tmax/(I**2) ) - 2*betaSquare )
	S_inv = S**(-1)
	return S_inv

# Parameters for air:
rho_air = 1.205e-3
I_air = 85.7e-6
Z_air = (6.+7.+8.+18.)/4.
A_air = (0.000124*12.011 + 0.755267*14.007 + 0.231781*15.999 + 0.012827*39.948)

#print Z_air, A_air

N=10000
Range_Al = np.zeros(20)
Range_con = np.zeros(20)
Range_air = np.zeros(20)
for k in range(20):
	Eincoming=np.linspace(0.2,10e3,20)
	E = np.linspace(0.1,Eincoming[k], N)
	InverseSP_Al = np.zeros(len(E))
	InverseSP_con = np.zeros(len(E))
	InverseSP_air = np.zeros(len(E))
	for i in range(N):
		InverseSP_Al[i] = InverseStoppingPower(E[i], rho_Al, I_Al, A_Al, Z_Al)
		InverseSP_con[i] = InverseStoppingPower(E[i], rho_con, I_con, A_con, Z_con)
		InverseSP_air[i] = InverseStoppingPower(E[i], rho_air, I_air, A_air, Z_air)
	Range_Al[k] = simps(InverseSP_Al, E)
	Range_con[k] = simps(InverseSP_con, E)
	Range_air[k] = simps(InverseSP_air, E)

plt.loglog(Eincoming, Range_Al)
plt.loglog(Eincoming, Range_con)
plt.loglog(Eincoming, Range_air)
plt.xlabel('Particle energy [MeV]')
plt.ylabel('Range [g/cm^2]')
plt.legend(['Aluminum', 'Portland concrete', 'Air'], loc=2)
plt.title('Particle range in different materials')
plt.grid(True, which="both")
plt.show()


#plt.show()



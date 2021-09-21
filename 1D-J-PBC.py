#!/usr/bin/env python
#use conda python 
import os
import numpy as np
from tqdm import tqdm 
import matplotlib.pyplot as plt
from matplotlib import patches as pt
rundir = os.getcwd()

#input parameters for Monte-carlo simulations
print("minimum temp set to 0.5 K; please set max temp in next line; 100 divisions are considered between min-max")
Temp_max = float(input("max. temp = "))
if Temp_max < 0.05:
   print("please enter a value greater than 0.06")

N = int(input("size of lattice = "))
J = float(input("exchange paramter = "))


#create lattice for the simulation
def hot_lattice(N):
    spin_value = np.array([-1, 1])
    lattice = np.random.choice(spin_value, size=N)
    return lattice

def cold_lattice(N):
    spin_value = np.array([1])
    lattice = np.random.choice(spin_value, size=N)
    return lattice

#calculate the magnetization
def Mag(lattice):
    Mag = np.sum(lattice)
    return Mag

#Calculate the energy of the lattice with PBC and no PBC
def calcEnergy_wPBC(lattice):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(lattice)):
            N = len(lattice)
            Si = lattice[i]
            Nj = lattice[(i+1)%N] + lattice[(i-1)%N]
            energy += -J*Si*Nj
    return energy/2.

def calcEnergy_woPBC(lattice):
    '''Energy of a given configuration'''
    energy = 0
    for i in range(len(lattice)):
            N = len(lattice)
            if i == 0:
               Si = lattice[i]
               Nj = lattice[(i+1)]
               energy += -J*Si*Nj
               print(Si*Nj)
            elif i == (N-1):
               Si = lattice[i]
               Nj = lattice[(i-1)]
               energy += -Si*Nj
               print(Si*Nj)
            elif i != (N-1) or i != 0.:
               Si = lattice[i]
               Nj = lattice[(i+1)%N] + lattice[(i-1)%N]
               energy += -Si*Nj
               print(Si*Nj)
    return energy/2.

# Monte carlo move
def mc(lattice, beta):
    for i in range(len(lattice)):
        N = len(lattice)
        r_n = np.random.randint(0,(N-1))
        Su = lattice[i]
        Nj = lattice[(i+1)%N] + lattice[(i-1)%N]
        dE = 2*J*Su*Nj
        if dE < 0.:
           Su *= -1
        elif np.random.uniform() < np.exp(-dE*beta):
           Su *= -1
        lattice[i] = Su
    return lattice

#temperature arrays
Temp_min = 0.05
T = np.linspace(Temp_min,Temp_max,30)
T_iter = int(len(T))

MC_step = 100
EQ_step = 100
 
#list to store final values
F_energys = []
F_Cps = []
F_Magnetizations = []
F_Chis = []
#Temperature = []

#Main code
for step in tqdm(range(T_iter)):
    lattice = cold_lattice(N)
    beta = 1.0/T[step]
    beta_sq = beta*beta
    for i in tqdm(range(EQ_step)):
        Energy_list = []
        Cp_list = []
        Magnetization_list = []
        Chi_list = []
         
        mc(lattice, beta)
        E = E_sq = M = M_sq = 0
        for i in range(MC_step):
            mc(lattice, beta)
            energy = calcEnergy_wPBC(lattice) 
            mag = abs(Mag(lattice))
            E += energy
            E_sq += energy**2
            M += mag
            M_sq += mag**2
        
        E_mean = E/MC_step
        E_sq_mean = E_sq/MC_step
        M_mean = M/MC_step
        M_sq_mean = M_sq/MC_step
        
        Energy = float(E_mean/float(N))
        Cp = beta**2 * ((E_sq_mean - E_mean**2)/float(N))
        Magnetization = M_mean/float(N)
        Chi = beta * ((M_sq_mean - M_mean**2)/float(N))
       
        #store the above values as list
        Energy_list.append(Energy); Cp_list.append(Cp); Magnetization_list.append(Magnetization); Chi_list.append(Chi)

    #to store list with temperature variation 
    F_energy = np.mean(Energy_list)
    F_energys.append(F_energy)
    F_Cp = np.mean(Cp_list)
    F_Cps.append(F_Cp)
    F_Magnetization = abs(np.mean(Magnetization_list))
    F_Magnetizations.append(F_Magnetization)
    F_Chi = np.mean(Chi_list)
    F_Chis.append(F_Chi)

f = plt.figure(figsize=(18, 10)); # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 );
plt.scatter(T, F_energys, s=30, marker='o', color='Red')
plt.xlabel("Temperature (T)", fontsize=10);
plt.ylabel("Energy", fontsize=10); 

sp =  f.add_subplot(2, 2, 2 );
plt.scatter(T, F_Magnetizations, s=30, marker='o', color='Blue')
plt.xlabel("Temperature (T)", fontsize=10); 
plt.ylabel("Magnetization ", fontsize=10);  

sp =  f.add_subplot(2, 2, 3 );
plt.scatter(T, F_Cps, s=30, marker='o', color='Green')
plt.xlabel("Temperature (T)", fontsize=10);  
plt.ylabel("Specific Heat ", fontsize=10); 

sp =  f.add_subplot(2, 2, 4 );
plt.scatter(T, F_Chis, s=30, marker='o', color='Cyan')
plt.xlabel("Temperature (T)", fontsize=10); 
plt.ylabel("Susceptibility", fontsize=10);  
      
  
plt.show()



# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
flex_code_string = """ 
TITLE 'Multiphysics Assignment 1'
COORDINATES cartesian1
VARIABLES
vy(threshold=0.1) !velocity in y
yd(threshold=0.1) !y-displacement
SELECT
ngrid = 1 
DEFINITIONS
!650000 = 800*(30*(1+m_fuel1/2000)+15*(1+m_fuel2/1000))+5*(m_fuel1)+5*(m_fuel2)
!m_fuel2=614000/17-m_fuel1 

m_fuel1=1*%s !Introductory fuel value
m_fuel2=614000/17-m_fuel1 !Sets initial fuel conditions, based on cost equation

qdot1=2000
vfuel1=2500 
m_engine1=320 
m_tank1=30*(1+m_fuel1/2000) !first engine and tank equations

qdot2=300
vfuel2=2800
m_engine2=150
m_tank2=15*(1+m_fuel2/1000) !second engine and tank equations

vi = 1 !initial velocity
m_payload=200 !payload of rocket

CD=0.2 !drag coefficient
A=1 !corresponding surface area
GC=6.674e-11 !gravitational constant
mEarth=5.9722e24 !mass of Earth
rEarth=6.3781e6 !Radius of Earth
Tnot=288.15 !Sea level standard temperature
pnot=101.325e3 !sea level standard atmospheric pressure
g=9.80665
M=0.0289644
R=8.31447
L=0.0065 !constants for fluidic drag force and forces

T1=Tnot-L*yd !temperature variations in atmosphere

cond=L*yd/Tnot !condition to check if it's in the atmosphere
p=if(cond>=1) then 0 else pnot*(1-cond)^(g*M/(R*L)) !accounts for leaving the atmosphere, rather than a negative rho

rhoAir=p*M/(R*T1) !when p is zero, out of the atmosphere.

tfuel1=m_fuel1/qdot1
tfuel2=m_fuel2/qdot2+tfuel1 !time it takes until fuel runs out

mfuel=if(t<tfuel1) then m_fuel1-qdot1*t+m_fuel2 else if (t<tfuel2) then m_fuel2-qdot2*(t-tfuel1) else 0 !fuel usage at each stage
m_tanks=if(t<tfuel1) then m_tank1+m_tank2 else if (t<tfuel2) then m_tank2 else 0 !assuming tank drops off after fuel is used up
m_engines=if(t<tfuel1) then m_engine1+m_engine2 else if (t<tfuel2) then m_engine2 else 0

mtotal=mfuel+m_engines+m_tanks+m_payload !total mass of each component

Ft1=qdot1*vfuel1 !thrust force of each component
Ft2=qdot2*vfuel2

Ft=if(t<tfuel1) then Ft1 else if (t<tfuel2) then Ft2 else 0 !changing thrust force

Fd=0.5*rhoAir*CD*A*(vy)^2 !changing drag force component

gforce=GC*mEarth*mtotal/(rEarth+yd)^2 !changing gravitational force for each mass

ay=(Ft-Fd-gforce)/mtotal !acceleration from thrust force

KE=1/2*mtotal*vy^2
PE=mtotal*g*(yd) !assuming potential is taken from the centre of the planet. 
TE=KE+PE !total energy to be kept in track of

INITIAL VALUES
vy = vi
yd = 0

EQUATIONS
vy: dt(vy) = ay
yd: dt(yd) = vy

BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START (0)   { Walk the domain boundary }
    LINE TO (1)
TIME 0 TO 100 halt (mfuel=0) and p=0 { if time dependent }

PLOTS
for t = 0 by endtime/100 to endtime
	!history(yd) at (0)
	history(yd) at (0) PrintOnly Export Format '#t#b#1' file ='Assign1FlexinOutput.txt'
    history(TE) at (0) PrintOnly Export Format '#t#b#1' file ='Assign1TEOutput.txt'
SUMMARY
report val (TE,0)
END"""

flex_filename = "tempflexfile.pde"

import subprocess
import numpy as np
import matplotlib.pyplot as plt

LowFuel=0
HighFuel=650000/17-LowFuel
NumPoints=18
Delta=HighFuel-LowFuel

BestTE=-1
BestFuel1=-1
BestFuel2=650000/17-BestFuel1

fuel_list= np.linspace(LowFuel,HighFuel,NumPoints)

for Fuel in fuel_list:
    with open(flex_filename,"w") as f:
        print(flex_code_string%Fuel,file=f) 
    subprocess.run(["C:\FlexPDE7\FlexPDE7.exe","-S",flex_filename],timeout=100)
    
    with open("tempflexfile_output\\Assign1FlexinOutput.txt") as f:
        flexoutputrawdata=np.loadtxt(f,skiprows=8)
    t = flexoutputrawdata[:,0]
    y = flexoutputrawdata[:,1]
    
   # print("With the fuel in tank one being %g, the rocket reached %g after %f seconds"%(Fuel,y[-1],t[-1]))
    
    with open("tempflexfile_output\\Assign1TEOutput.txt") as f:
        flexoutputrawdata=np.loadtxt(f,skiprows=8)
    t = flexoutputrawdata[:,0]
    TE = flexoutputrawdata[:,1]
    TEfinal=TE[-1]
    print("With the fuel in tank one being %g, the rocket reaches a total energy of %g after %f seconds"%(Fuel,TE[-1],t[-1]))
    plt.plot(t,TE)
    
    if(TEfinal>BestTE):
        BestTE=TEfinal
        BestFuel1=Fuel
        BestFuel2=650000/17-BestFuel1
    
#plt.legend(fuel_list)
#plt.show()
    
    
        
print(f"Best fuel1 amount was {BestFuel1} with a final TE of {BestTE} J")
    
    
# -*- coding: utf-8 -*-


# Version 1- basic modelling-final spectral index plot only
#Ver 1.0.1 added diagnostics plots and tables when modelling is finished
#Ver 1.0.2 adding diagnostic plots to every 10000 iterations so program doesn't need to be run several times
#Ver 1.1.0 version that is used in Mulcahy et al, 2016
#Ver 1.1.1 version has a user parset interface

# import modules and functions
from __future__ import division
import numpy as np
from diffunction import *
import os
import sys
import csv
from optparse import OptionParser


VERSION = '1.1.1'

####proper plots and output textfiles
#### This code models the diffusion loss equation


#############################################PARSER OPTIONS#########################################################

if __name__ == '__main__':
    
    """ Handle all parsing here if started from the command line"""
    parser = OptionParser(usage="%prog <input parameter file>",  version="%prog " + VERSION)


    (options, args) = parser.parse_args()
    
    print "CREPE ver. " + VERSION
    print "Written by David D. Mulcahy"
    print ""
    print "Parsing parameter file..."
    params = readfromparset(args[0])


D = params.diffcoeffr
Dz = params.diffcoeffz
kappa = 1
nt = params.iter
outputarray = 'output.npy'


####Creating plot directories####

if os.path.exists(params.outputdir):
    os.system('rm -rf '+str(params.outputdir))
    os.system('mkdir '+str(params.outputdir))
else:
    os.system('mkdir '+str(params.outputdir))


#Parameters to be used

####SPATIAL #####
startr = 0.05
endr = 15.05
dr = 0.05
nr = (endr - startr)/dr
r = np.linspace(startr,endr,nr+1)*3.08567758e21 #Initialize spatial axis
r0 = 4.64394476e22 #10.05kpc in cms
rprime = r/r0 #scaled radius
drprime = rprime[1]-rprime[0]


# creating r for derivative only
derivstartr = -0.05
derivendr = 15.15
derivenr = (derivendr - derivstartr)/dr
deriver = np.linspace(derivstartr,derivendr,derivenr+1)*3.08567758e21 #Initialize spatial axis
deriverprime = deriver/r0

##### TIME #####
dt = 0.01 # timestep 
t0=3.1556926e13 # million years in seconds 

### ENERGY ###
starte=0.01 #in GeV
ende=5
ne=300
E = np.linspace(starte,ende,ne+1)

#Energy Dependence of diffusion
#n=0
#Diff = np.ones(len(E))
#ENERGY=E/3
#print ENERGY
#for finddiff in ENERGY:
#    if ENERGY[n]<=1:
#        Diff[n]=D
#        n=n+1
#    else:
#        Diff[n]=D*(ENERGY[n])**kappa
#        n=n+1


#print'DIFF is', Diff

           
### Coefficients for injection spectrum #####
K = 1
gamma_0 = -2.001 #injection energy slope -2 is equal to -0.5 in spectral index

#######create Source Term#############
N = K*(E**(gamma_0)) ###!!!!N.B linear E must go into this function!!!!#####

############Input RADIAL Magnetic Field####################
bfeldy = np.loadtxt(params.inputmagnetic,unpack=True, usecols=[1])
print bfeldy
bfeld0 = np.max(bfeldy) #normalise the magnetic field
print 'Bfeld0 is',bfeld0 
bfeldprime = bfeldy/bfeld0 # making bfield dimensionless

########Input radial SNR for injection##################
sourcey = np.loadtxt(params.inputsource,unpack=True, usecols=[1])

inputQ = sourcey*0.02
Q0 = np.max(inputQ)/2
inputQ = (inputQ/Q0)

# checking if the second derivative at nr is 0, all points in this area should be zero
if inputQ[nr]==0 and inputQ[nr-1]==0 and inputQ[nr-2]==0:
    print 'initial source boundary condition is okay'
else:
    print 'initial source boundary condition unacceptable, the final three entries in input source function need to be zero. Please modify your source function.'
    sys.exit('Terminating program')
    
n=0
Q = np.zeros((301,301))
for work in N:
  Q[n,:]=N[n]*inputQ  #fullarray[n,:] is radial, fullarray[:,n] is energy
  n=n+1

print Q  

###############Creating coefficients for diffusion and energy losses#############
E0 = np.max(E)
E0 = 1
Eprime = E/E0
dEprime = np.diff(Eprime)
print 'E0 is', E0
print 'E0 in ergs is', E0*1.6e-12
print 't0 is',t0
print 'r0 is', r0
print 'D is', D
#theta = (t0*Diff)/(r0**2) #turn on energy dependent diffusion
theta = (t0*D)/(r0**2)
print 'Theta is',theta


alpha = 8e-17
xi = 2e-6
phi = alpha*t0*E0
sigma = xi*t0*E0*(bfeld0)**2
print 'phi is',phi
print 'sigma is',sigma

diffusioncourant = (theta*dt)/drprime
print 'Diffusion courant condition is', diffusioncourant


energysigmacourant = (sigma*dt)/dEprime
print 'Sigma courant condition is', energysigmacourant
print dEprime


#if np.all(diffusioncourant) < 1.0:
#    print 'Diffusion Courant Condition has passed'
#else:
#    print 'Diffusion Courant Condition has not passed'
#    sys.exit('Terminating program')


#Computing the escape times in radius      
#scaleheight[:] = 1.4 #uncomment if one wants Mao et al 2015 scale heights
scaleheight = scaleheight*4 #magnetic field scale height is 4 times the synchrotron scale height
scaleheightdif = ((scaleheight*3.08567758e21)**2)/(Dz) #convert to cms and find the escape time independent of energy


#######RADIAL DIFFUSION SECTION---PLEASE CHECK HERE FOR BOUNDARY CONDITIONS#######
#Spatial part of Right Hand Side of Equation in cylindrical coordinates
def rhs(u):
  n = 0
  newQE = np.zeros((301,301))
  for cyclefinite in u[n,:]:
    Dhigh = u[n, nr] + (u[n, nr]-u[n, nr-1])
    D2high = u[n, nr] + (u[n, nr]-u[n, nr-2])
    Dlow = D2low = u[n,0]
    uplus = np.hstack((Dlow,D2low,u[n,:],Dhigh,D2high))
    uplusdev= (-1*np.roll(uplus,-2) + 16*np.roll(uplus,-1) - 30*(uplus) + 16*np.roll(uplus,1) - np.roll(uplus,2))/( 12*(drprime**2) ) # centered differences on all points
    newQE[n,:] =uplusdev[2:-2] #get rid of the extrapolated points
    n=n+1
  return newQE


# Energy part of the Right Hand Side of Equation 
def rhsenergy(u):
	n=0		
	newQE = np.zeros((301,301))
	for cyclefiniteenergy in u[:,n]:
		QElow=findinterpolpointslowenergy(u[:,n]) 
		Elow=findinterpolpointslowenergy(Eprime)
		QEhigh=findinterpolpointshighenergy(u[:,n],ne)
		Ehigh=findinterpolpointshighenergy(Eprime,ne)
		newQE[1:,n] = ((((np.roll(Eprime,-1)**2)*np.roll(u[:,n],-1)) - ((Eprime**2)*u[:,n]))[1:]/dEprime) 
		newQE[0,n] = (((Eprime[0]**2)*u[0,n]) - ((Elow**2)*QElow))/(Eprime[0]- Elow) 
		newQE[ne,n] = (((Ehigh**2)*QEhigh) - ((Eprime[ne]**2)*u[ne,n]))/(Ehigh - Eprime[ne])
		#if n==200:
		#	print newQE[ne,n]
		n=n+1  
	return newQE


#Adding all three parts of right hand side
def fullrightside(u):
  finalu = K*Q + theta*rhs(u) + bfeldprime*sigma*rhsenergy(u) + phi*rhsenergy(u)
  return finalu

if params.inputarray:
    print 'Input numpy array accepted'
    u = np.load(params.inputarray)
else:
    print 'No Input inserted'
    print 'Will begin modelling naturally'
    u=K*Q #start with a single injection

i=1

#factors to scale back radius and energy
finalr=(rprime*r0)/3.08567758e21 # real radius in kpc for plots
finale=(Eprime*E0)/1.6e-12 # real energy in ev for plots
finale=finale/1e09

#initialise arrays for diagnostics 
timearray = np.array([])
dianosticcenter = np.array([])
dianosticinter = np.array([])
dianosticarm = np.array([])
dianosticextended = np.array([]) 
allboundarypoints = np.array([])
time = np.array([])
totalcrefile = open('creevolution.txt','w') # open up a file to show the total number of CRe in model at each iteration

# Runga-Kutta Loop

while i<nt:
	k1 = dt*fullrightside(u)
	k2 = dt*fullrightside(u+k1/2.0)
	k3 = dt*fullrightside(u+k2/2.0)
	k4 = dt*fullrightside(u+k3)
	u = u+(k1+2*k2+2*k3+k4)/6.0
	u = u.clip(min=0)
        if i%1000==0:
            print i
            #print u[0,:]
            #print 'creating prelimnary plots'
            finaloutput50MHZ, LOFARLBAenergy = find_nearest(u,5.0e+01,bfeldy,E)
            finaloutput151MHZ, LOFARenergy = find_nearest(u,1.51e+02,bfeldy,E)
            finaloutput330MHZ, GMRT330energy = find_nearest(u,3.3e+02,bfeldy,E)
            finaloutput610MHZ, GMRT610energy = find_nearest(u,6.1e+02,bfeldy,E)
            finaloutput1400MHZ, VLAenergy = find_nearest(u,1.4e+03,bfeldy,E)
            finaloutput2600MHZ, SBANDenergy = find_nearest(u,3.0e+03,bfeldy,E)
            modelspecHBAVLA = (np.log10(finaloutput1400MHZ)-np.log10(finaloutput151MHZ))/(np.log10(VLAenergy)-np.log10(LOFARenergy))
            freqspecHBAVLA = (-1*((-1*modelspecHBAVLA)-1))/2
            modelspecHBALBA = (np.log10(finaloutput151MHZ)-np.log10(finaloutput50MHZ))/(np.log10(LOFARenergy)-np.log10(LOFARLBAenergy))
            freqspecHBALBA = (-1*((-1*modelspecHBALBA)-1))/2
            createfinalplot(finalr,freqspecHBAVLA,i)
            createfinalplotHBALBA(finalr,freqspecHBALBA,i)
            outputfiles(finalr,freqspecHBAVLA,freqspecHBALBA,finaloutput50MHZ,finaloutput151MHZ,finaloutput330MHZ,finaloutput610MHZ,finaloutput1400MHZ,finaloutput2600MHZ,i)
            totalN=np.sum(u)
            totalcrefile.write(str(i) + " " + str(totalN) + "\n")
    	i=i+1

print timearray.shape
print dianosticarm.shape
print 'simulation finished'
print 'creating outputs and fits'
print u
np.save(str(outputarray)+'.npy',u)

#creating final plots for analysis
finaloutput50MHZ, LOFARLBAenergy = find_nearest(u,5.0e+01,bfeldy,E)
finaloutput151MHZ, LOFARenergy = find_nearest(u,1.51e+02,bfeldy,E)
finaloutput330MHZ, GMRT330energy = find_nearest(u,3.3e+02,bfeldy,E)
finaloutput610MHZ, GMRT610energy = find_nearest(u,6.1e+02,bfeldy,E)
finaloutput1400MHZ, VLAenergy = find_nearest(u,1.4e+03,bfeldy,E)
finaloutput2600MHZ, SBANDenergy = find_nearest(u,3.0e+03,bfeldy,E)
modelspecHBAVLA = (np.log10(finaloutput1400MHZ)-np.log10(finaloutput151MHZ))/(np.log10(VLAenergy)-np.log10(LOFARenergy))
freqspecHBAVLA = (-1*((-1*modelspecHBAVLA)-1))/2
modelspecHBALBA = (np.log10(finaloutput151MHZ)-np.log10(finaloutput50MHZ))/(np.log10(LOFARenergy)-np.log10(LOFARLBAenergy))
freqspecHBALBA = (-1*((-1*modelspecHBALBA)-1))/2
createfinalplot(finalr,freqspecHBAVLA,i)
createfinalplotHBALBA(finalr,freqspecHBALBA,i)
outputfiles(finalr,freqspecHBAVLA,freqspecHBALBA,finaloutput50MHZ,finaloutput151MHZ,finaloutput330MHZ,finaloutput610MHZ,finaloutput1400MHZ,finaloutput2600MHZ,i)
totalN=np.sum(u)
totalcrefile.write(str(i) + " " + str(totalN) + "\n")

print 'Finished after '+str(i)+' iterations'
sys.exit()

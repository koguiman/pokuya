###This script is used to generate the pokuya software
###pokuya is an open source software developed to obtain 
###the same basic functionality as a power quiality analyzer

import pylab
import numpy, scipy
import matplotlib.pyplot as plt 
from math import *
from cmath import *
from scipy import linalg as nn
import numpy as np
import csv
from scipy import signal as sg
from csvtonum import *

##########Single phase tests
###########inputs part
name='Test3.csv' 
fe=float(50) #electrical frequecy. Eur=50Hz, American countries=60Hz
Ncycle=3 #number of cycles to be analyzed
###########output ports
reportname='{0}{1}'.format('reportof_',name)
########################
##open the file with data
x=csv.reader(open(name),delimiter=',')
###load data and info from source file
[data,info]=csvtodato(x)
####read data
time=data[:,0]; #time signal from scopes or pcs
time=time-min(time) #time signal shifted to zero t
Ts=abs(time[1]-time[0]); #sampling period
I=data[:,1] #current signal
V=data[:,3] #voltage signal

###############If a resampling is necessary

Fs=ceil(1.0/(Ts)); #sampling frequecy
Ln=I.shape[0] #size of the signals
df=(Fs/Ln); ##steps for the frequency domain, dato from the normal signals
###Remark: in case df is not multiple of fe=50 Hz it is necessay the resampling
if (ceil(fe/df)!=fe/df):
        Ttotal2=Ncycle*(1./fe);
        df2=1./Ttotal2
        l1=ceil(Ttotal2/Ts)
        msamp=100
        Ts2=(1./fe/msamp)
        Ln2temp=ceil(Ttotal2/Ts2)
        Ln2=Ln2temp
        I,time=sg.resample(I[0:l1], Ln2, time[0:l1], axis=0, window=None)
        V,time=sg.resample(V[0:l1], Ln2, time[0:l1], axis=0, window=None)
        Ln=Ln2
        Fs=1./Ts2#new sampling frequecy
###########Fourier transforms part
Idfft=2*np.fft.fft(array(I))/Ln
Vdfft=2*np.fft.fft(array(V))/Ln
##Fs=1.0/(Ts); #sampling frequecy
Freq = (Fs/Ln)*linspace(0,Ln-1,num=Ln);

F_Id=abs(Idfft) #Fourier magnitude Current
ph_Id=np.arctan2((Idfft).real,(Idfft).imag); #angle 
F_Vd=abs(Vdfft) ##Fourier  magnitude voltage
ph_Vd=np.arctan2((Vdfft).real,(Vdfft).imag); #angle  voltage

#################################
pinst=V*I #instantaneous power
################################
###########Report stage
I0=F_Id[0] #DC component
V0=F_Vd[0] #DC component
f_50=ceil(fe/df2); #fundamental frequency position
I50=F_Id[f_50]/sqrt(2.)#fundamental current component rms 
V50=F_Vd[f_50]/sqrt(2.)#fundamental voltage component
maxhar=50; #maximmun harmonic to be represented
Ihs=F_Id[f_50:f_50*maxhar:f_50]/sqrt(2.)#Harmonics current rms
Vhs=F_Vd[f_50:f_50*maxhar:f_50]/sqrt(2.)#Harmonics current rms
Fhs=Freq[f_50:f_50*maxhar:f_50]#Harmonics value 
I_THDf=sqrt(sum(Ihs[1:Ihs.shape[0]]**2))/I50 #THD respect the fundamental
V_THDf=sqrt(sum(Vhs[1:Vhs.shape[0]]**2))/V50 #THD respect the fundamental
I_THDr=I_THDf/sqrt(1+(I_THDf**2))# THD respect the distorted signal
V_THDr=V_THDf/sqrt(1+(V_THDf**2))# THD respect the distorted signal
DPF=1./sqrt(1.+I_THDf) #displacement power factor 
Irms=I50*sqrt(1+I_THDf**2)
Vrms=V50*sqrt(1+V_THDf**2)
PF=DPF*(I50/Irms) ##power factor
cos_phi=cos(ph_Vd[f_50]-ph_Id[f_50]) #cos? phi, phi:angle between the voltage and current at fe
Paver=mean(pinst) #average power
##put more analysis
##################################################################
###Report file generation
with open(reportname, 'wb') as csvfile:
    stwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    stwriter.writerow(['Report of :'] + [name])
    stwriter.writerow(['Fundamental frequency [Hz]'] + [str(fe)])
    stwriter.writerow(['DC current [A]'] + [str(I0)])
    stwriter.writerow(['DC voltage [V]'] + [str(V0)])
    stwriter.writerow(['Fundamental current rms [A]'] + [str(I50)]+['Fundamental current peak [A]'] + [str(I50*sqrt(2))])
    stwriter.writerow(['Fundamental voltage rms [V]'] + [str(V50)]+['Fundamental voltage peak [V]'] + [str(V50*sqrt(2))])
    stwriter.writerow(['Current: THDf'] + [str(I_THDf)])
    stwriter.writerow(['Voltage: THDf'] + [str(V_THDf)])
    stwriter.writerow(['Current: THDr'] + [str(I_THDr)])
    stwriter.writerow(['Voltage: THDr'] + [str(V_THDr)])
    stwriter.writerow(['Displacement power factor DPF'] + [str(DPF)])
    stwriter.writerow(['Power factor PF'] + [str(PF)])
    stwriter.writerow(['Harmonics']+['Ih']+['%I1']+['Vh']+['%V1'])
    # Fhs,Ihs,Vhs
    Harmo=c_[Fhs,Ihs,(100*Ihs/I50),Vhs,(100*Vhs/V50)]
    for i in Harmo:
            stwriter.writerow([str(i)])
    
#################################################################
###########Ploting part

####time domain signals
plt.figure(1)
plt.subplot(311)
plt.plot(time[0:Ln],I)
plt.title('Current [A]')
plt.ylabel('I(t) [A]')
plt.grid(True)
plt.subplot(312)
plt.plot(time[0:Ln],V)
plt.title('Voltage [V]')
plt.ylabel('V(t) [V]')
plt.grid(True)
#plt.xlabel('Time [s]')
plt.subplot(313)
plt.plot(time[0:Ln],pinst)
plt.title('Instantaneous power [W]')
plt.ylabel('P(t) [W]')
plt.xlabel('Time [s]')
plt.grid(True)

###Fourier spectrum
lim=int(Freq.shape[0]/2)
plt.figure(2)
plt.subplot(211)
#plt.plot(Freq[0:lim],F_Id[0:lim])
plt.bar(np.hstack([array([0.0]),Fhs]),np.hstack([array([I0]),Ihs]), width=0.84, label="Current",align="center")
plt.xlim([0,50*50])
plt.ylabel('I($\omega$) [A$_{rms}$]')
plt.grid(True)
plt.subplot(212)
#plt.plot(Freq[0:lim],F_Vd[0:lim])
plt.bar(np.hstack([array([0.0]),Fhs]),np.hstack([array([V0]),Vhs]), width=0.84, label="Current",align="center")
plt.xlim([0,50*50])
plt.ylabel('V($\omega$) [V$_{rms}$]')
plt.grid(True)
###Fourier spectrum in % of the I50
plt.figure(3)
plt.subplot(211)
#plt.plot(Freq[0:lim],F_Id[0:lim])
plt.bar(np.hstack([array([0.0]),Fhs]),100*np.hstack([array([I0]),Ihs])/I50, width=0.84, label="Current",align="center")
plt.xlim([0,50*50])
plt.ylim([0,101])
plt.ylabel('I($\omega$) [%]')
plt.grid(True)
plt.subplot(212)
#plt.plot(Freq[0:lim],F_Vd[0:lim])
plt.bar(np.hstack([array([0.0]),Fhs]),100*np.hstack([array([V0]),Vhs])/V50, width=0.84, label="Current",align="center")
plt.xlim([0,50*50])
plt.ylim([0,101])
plt.ylabel('V($\omega$) [ %]')
plt.grid(True)
plt.xlabel('Frequency [Hz]')

plt.show()

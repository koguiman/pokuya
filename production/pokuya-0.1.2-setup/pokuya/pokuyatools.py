import pylab,scipy,csv,cmath,math
from math import *
from numpy import *
from scipy import signal as sg
import numpy as np
import matplotlib.pyplot as plt


#singlephase analysis
def sphan(X):
    ###  X:is a dictionary with the time, voltage and current signals for a single phase analysis
    data=X['data'];    Ncycle=X['Ncycle']; fe=X['fe']; channcurr=X['channcur']
    time=data[:,0]; #time signal from scopes or pcs
    time=time-min(time) #time signal shifted to zero t
    Ts=abs(time[1]-time[0]); #sampling period
    if (channcurr==1):
        I=data[:,1] #current signal
        V=data[:,3] #voltage signal
    elif (channcurr==2):
        I=data[:,3] #current signal
        V=data[:,1] #voltage signal
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
    Iph_hs=ph_Id[f_50:f_50*maxhar:f_50]*(180./np.pi)#Harmonics current angle in deg
    Vhs=F_Vd[f_50:f_50*maxhar:f_50]/sqrt(2.)#Harmonics voltage rms
    Vph_hs=ph_Vd[f_50:f_50*maxhar:f_50]*(180./np.pi)#Harmonics voltage angle in deg
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
    Harmo=c_[Fhs,Ihs,(100*Ihs/I50),Iph_hs,Vhs,(100*Vhs/V50),Vph_hs]
    colinf=array(['Freq Hz','I_harmonics rms','I_harmonics %','angle_I deg','V_harmonics rms','V_harmonics %','angle_V deg'])
    Harmoinf={'Harmo':Harmo,'colinfo':colinf}
    current={'I0':I0,'I50':I50,'Irms':Irms,'I_THDf':I_THDf,'I_THDr':I_THDr}
    voltage={'V0':V0,'V50':V50,'Vrms':Vrms,'V_THDf':V_THDf,'V_THDr':V_THDr}
    ##put more analysis
    report={'Fhs':Fhs,'current':current, 'voltage':voltage,'fe':f_50,'Harmonics':Harmoinf,'DPF':DPF,'PF':PF,'cos_phi':cos_phi,'pinst':pinst,'time':time,'I':I,'V':V}
    return report


def savereport(report,reportname,name):
    with open(reportname, 'wb') as csvfile:
        stwriter = csv.writer(csvfile, delimiter=' ',  quoting=csv.QUOTE_MINIMAL)#quotechar= '|',
        stwriter.writerow(['Report of :'] + [name])
        stwriter.writerow(['Fundamental frequency [Hz]'] + [str(report['fe'])])
        stwriter.writerow(['DC current [A]'] + [str(report['current']['I0'])])
        stwriter.writerow(['DC voltage [V]'] + [str(report['voltage']['V0'])])
        stwriter.writerow(['Fundamental current rms [A]'] + [str(report['current']['I50'])]+['Fundamental current peak [A]'] + [str(report['current']['I50']*sqrt(2))])
        stwriter.writerow(['Fundamental voltage rms [V]'] + [str(report['voltage']['V50'])]+['Fundamental voltage peak [V]'] + [str(report['voltage']['V50']*sqrt(2))])
        stwriter.writerow(['Current: THDf'] + [str(report['current']['I_THDf'])])
        stwriter.writerow(['Voltage: THDf'] + [str(report['voltage']['V_THDf'])])
        stwriter.writerow(['Current: THDr'] + [str(report['current']['I_THDr'])])
        stwriter.writerow(['Voltage: THDr'] + [str(report['voltage']['V_THDr'])])
        stwriter.writerow(['Displacement power factor DPF'] + [str(report['DPF'])])
        stwriter.writerow(['Power factor PF'] + [str(report['PF'])])
        stwriter.writerow(report['Harmonics']['colinfo'])
        # Fhs,Ihs,Vhs
        for i in report['Harmonics']['Harmo']:
                stwriter.writerow([str(i)])
    
#################################################################
def grafic(report):
    time=report['time']
    I=report['I']
    V=report['V']
    pinst=report['pinst']
    Freq=report['Fhs']
    Ln=I.shape[0] #size of the signals
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
    plt.figure(2)
    plt.subplot(211)
    #plt.plot(Freq[0:lim],F_Id[0:lim])
    plt.bar(np.hstack([array([0.0]),report['Fhs']]),np.hstack([array([report['current']['I0']]),report['Harmonics']['Harmo'][:,1]]), width=0.84, label="Current",align="center")
    plt.xlim([0,50*50])
    plt.ylabel('I($\omega$) [A$_{rms}$]')
    plt.grid(True)
    plt.subplot(212)
    #plt.plot(Freq[0:lim],F_Vd[0:lim])
    plt.bar(np.hstack([array([0.0]),report['Fhs']]),np.hstack([array([report['voltage']['V0']]),report['Harmonics']['Harmo'][:,4]]), width=0.84, label="Current",align="center")
    plt.xlim([0,50*50])
    plt.ylabel('V($\omega$) [V$_{rms}$]')
    plt.grid(True)
    ###Fourier spectrum in % of the I50
    plt.figure(3)
    plt.subplot(211)
    #plt.plot(Freq[0:lim],F_Id[0:lim])
    plt.bar(np.hstack([array([0.0]),report['Fhs']]),100*np.hstack([array([report['current']['I0']]),report['Harmonics']['Harmo'][:,1]])/report['current']['I50'], width=0.84, label="Current",align="center")
    plt.xlim([0,50*50])
    plt.ylim([0,101])
    plt.ylabel('I($\omega$) [%]')
    plt.grid(True)
    plt.subplot(212)
    #plt.plot(Freq[0:lim],F_Vd[0:lim])
    plt.bar(np.hstack([array([0.0]),report['Fhs']]),100*np.hstack([array([report['voltage']['V0']]),report['Harmonics']['Harmo'][:,4]])/report['voltage']['V50'], width=0.84, label="Voltage",align="center")
    plt.xlim([0,50*50])
    plt.ylim([0,101])
    plt.ylabel('V($\omega$) [ %]')
    plt.grid(True)
    plt.xlabel('Frequency [Hz]')

    plt.show()

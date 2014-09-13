import scipy,numpy
from scipy import signal as sgn

#singlephase analysis
def sphan(X):
    ###  X:is a dictionary with the time, voltage and current signals for a single phase analysis
    data=X['data'];    Ncycle=X['Ncycle']; fe=X[fe]; chnncurr=X['chancur']
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
    Harmo=c_[Fhs,Ihs,(100*Ihs/I50),Vhs,(100*Vhs/V50)]
    colinf=array(['Freq','I_harmonics rms','I_harmonics %','V_harmonics rms','V_harmonics %'])
    Harmoinf={'Harmo':Harmo,'colinfo':colinf}
    current={'I0':I0,'I50':I50,'Irms':Irms,'I_THDf':I_THDf,'I_THDr':I_THDr}
    voltage={'I0':V0,'I50':V50,'Irms':Vrms,'I_THDf':I_THDf,'V_THDr':V_THDr}
    ##put more analysis
    report={'current':current, 'Voltage':voltage,'f_50':f_50,'Harmonics':Harmoinf,'DPF':DPF,'PF':PF,'cos_phi':cos_phi}
    return report


def savereport(report,reportname):
    with open(reportname, 'wb') as csvfile:
    stwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    stwriter.writerow(['Report of :'] + [name])
    stwriter.writerow(['Fundamental frequency [Hz]'] + [str(fe)])
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
    for i in report['Harmonics']'Harmo']:
            stwriter.writerow([str(i)])
    
#################################################################
    print 'nothing'

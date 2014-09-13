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
from pokuyatools import *

##########Single phase tests
###########inputs part
name='Test3.csv' 
fe=float(50) #electrical frequecy. Eur=50Hz, American countries=60Hz
Ncycle=3 #number of cycles to be analyzed
chnncurr=1
###########output ports
reportname='{0}{1}'.format('reportof_',name)
########################
##open the file with data
x=csv.reader(open(name),delimiter=',')
###load data and info from source file
[data,info]=csvtodato(x)

#create the object X for analysis
X={'data':data,'chancur':chnncurr,'Ncycles':Ncycle,'fe':fe}
############################
report=sphan(X)

##################################################################
###Report file generation
savereport(report,reportname)
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

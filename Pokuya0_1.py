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
from formlayout import fedit

from csvtonum import *
from pokuyatools import *

###GUI
###########################
datalist = [("Provide file name",""),
            ('Number of cycles',3),
            ('Save report',True),
            ('Current channel',1)
            ]

filename=fedit(datalist,title="Pokuya",comment="Please provide")


##########Single phase tests
###########inputs part
name=str(filename[0])
fe=float(50) #electrical frequecy. Eur=50Hz, American countries=60Hz
Ncycle=filename[1] #number of cycles to be analyzed
saverep=filename[2]
channcur=filename[3]


###########output ports
reportname='{0}{1}'.format('reportof_',name)
########################
##open the file with data
x=csv.reader(open(name),delimiter=',')
###load data and info from source file
[data,info]=csvtodato(x)

#create the object X for analysis
X={'data':data,'channcur':channcur,'Ncycle':Ncycle,'fe':fe}
############################
report=sphan(X)

##################################################################
###Report file generation
if (saverep==True):
    savereport(report,reportname,name)

#################################################################
###########Ploting part
grafic(report)

import csv
from numpy import *
from scipy import *
import numpy as np

name='Test1.csv'
x=csv.reader(open(name),delimiter=',')
ind = 0
a=0
data=zeros((1,4))
info=[];
for row in x:
    ind+=1
    if ind<18:
          info.insert(ind,','.join(row));
    if ind>18:
        b=double( np.hstack((row[3:5],row[9:11])))
        if a>0 and b.shape[0]>0:
            data=vstack((data,b))
        elif a==0:
            data[a,:]=b
        a+=1


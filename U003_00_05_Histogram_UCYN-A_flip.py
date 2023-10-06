'''
Created on Sep 30, 2021

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from sf import *

a = genfromtxt('../data/N2fixData.csv',delimiter=',')
UCYNA = genfromtxt('../data/UCYN-Adata.csv',delimiter=',')

figure(1)
hist([[0],[0]],edgecolor='k',color=['pink','green'],label=['Others','UCYN-A'])
hist([UCYNA,a],edgecolor='k',color=['green','pink'],stacked=True)
xlabel('N$_{2}$ fixation rate (mol N mol N$^{-1}$d$^{-1}$)')
ylabel('Frequency')
#hist(a,bins=arange(0,13,0.5))

legend(edgecolor='k')

sf('','histogram',450)

show()
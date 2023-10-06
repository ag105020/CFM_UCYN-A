'''
Created on Sep 30, 2021

@author: keiin
'''

from pylab import *
from FigSetting2 import *
from sf import *
from E_BiosynthesisFromNH4 import *

Taniguchi = genfromtxt('../data/Taniguchi.csv',delimiter=',').T
Maranon = genfromtxt('../data/Maranon_full.csv',delimiter=',').T
rcParams['figure.figsize']=3.5,2.844
rcParams.update({'font.size': 9,
                 'lines.markersize':8,
                 'lines.markeredgewidth':1})
rcParams.update({'ytick.major.pad': 1})
rcParams.update({'xtick.major.pad': 1})
rcParams.update({'axes.linewidth':1}) 

ep = 0.6
E = E_BiosynthesisFromNH4(ep)

figure(1)
plot(10**Taniguchi[0],(Taniguchi[1]*(1+E)),'o')
xlabel('V ($\mu$m$^{-3}$)')
ylabel('P$_{C}^{max}$ (d$^{-1}$)')
xscale('log')
xticks([10,100])
ylim(bottom=0,top=2.8)
#ylim(bottom=0,top=1.6+1e-10)
#yticks(arange(0,1.6+1e-1,0.2))
sf('','Taniguchi',450)

figure(2)
plot(Maranon[1],Maranon[2]*12,'o')
xlabel('V ($\mu$m$^{-3}$)')
ylabel('P$_{C}^{max}$ (d$^{-1}$)')
xscale('log')
print(nanmax(Maranon[2]*12))
ylim(bottom=0,top=2.8)
a = (4.29,1460.6)
b = array([0,0])
c = array([3,3])
print(size(c >= b))

fill_between(a,b,c, where=c >= b, facecolor='#e6e6e6',edgecolor = "none")
sf('','Maranon',450)

show()
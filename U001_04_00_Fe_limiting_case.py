'''
Created on Feb 17, 2020

@author: keiin
'''



from pylab import *
from E_BiosynthesisFromNH4 import *
from N2fixBudget import *
from FigSetting2 import *
from sf import *

A = 0.216    #(various dimensions bassed on B) Factor converting V to Qc
B = 0.939    #(dimensionless) Power factor (Menden-Deuer and Lessard 2000)
Ycn = 6.6    #(molC molN-1) C:N
E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate
Apho = 1.2*(1+E)     #(d-1) Photosythesis rate divided by cellular C
YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
#YresN2fix = 0
Y = YcnN2fix+YresN2fix
C = (1+E+Y)       # Case without energy saving by light harvesting by UCYN-A
#C = (1+E+Y/2*1.5) # Case with energy saving by light harvesting by UCYN-A

An2fix = 5.5*0.87*2.3272749437058176   #value based on DDA *5.5 times their own. 2.32... is the (Richelia N/heterocyst N) Value from D002 03 00

An2fix3 = 15.5

FF = 0.3  #Fe liimtation factor

An2fix4 = An2fix3*FF
Apho4 = Apho*FF

RhRd = arange(0.1,10+0.1,0.1)   #(um/um) Radius of haptophyto / radius of diazotroph

MuC = Apho*(RhRd)**(3*B)/(C*(RhRd**(3*B)+1))
MuC4 = Apho4*(RhRd)**(3*B)/(C*(RhRd**(3*B)+1))

MuN = An2fix/(1+RhRd)**(3*B)
MuN3 = An2fix3/(1+RhRd)**(3*B)
MuN4 = An2fix4/(1+RhRd)**(3*B)

Dpi = 300

figure(1)
plot((2.5,2.5),(0,0.7),':k')
plot(RhRd,MuC,'b',label='C lim: Replete')
plot(RhRd,MuN3,'r',label='N lim: Replete')
plot(RhRd,MuC4,'--b',label='C lim: Fe lim')
plot(RhRd,MuN4,'--r',label='N lim: Fe lim')
ylim(0,0.7)
xlim(left=0)
xlabel('R$_{H}$/R$_{D}$')
ylabel('$\mathit{\mu}$ (d$^{-1}$)')
legend(fontsize=17,loc=4,bbox_to_anchor=(1,0.3))
sf('','Felim',Dpi)

show()
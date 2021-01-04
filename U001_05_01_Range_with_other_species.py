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
Ycn = 6.3    #(molC molN-1) C:N
E = E_BiosynthesisFromNH4(0.6) #(mol C mol C-1) CO2 production rate
Apho = 1.26*(1+E)     #(d-1) Photosythesis rate divided by cellular C
Apho2 = 0.76*(1+E)
Apho3 = 1.01*(1+E)

YcnN2fix,YresN2fix = N2fixBudget(0.6) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
#YresN2fix = 0
Y = YcnN2fix+YresN2fix
C = (1+E+Y)       # Case without energy saving by light harvesting by UCYN-A
C = (1+E+Y/2*1.5) # Case with energy saving by light harvesting by UCYN-A
#An2fix = 6.1/60*12 #(d-1) N2 fixation rate per N (based on 6.1 fmol h-1 cell-1 and fmol h-1 cell-1)
An2fix = 5*(18.845*4+61.0125)/61.0125   #value based on DDA
An2fix = 5.5*0.87*2.3272749437058176   #value based on DDA *5.5 times their own. 2.32... is the (Richelia N/heterocyst N) Value from D002 03 00

print(An2fix)

An2fix2 = 0.8  #for Crocosphaera 
An2fix3 = 15.5
An2fix4 = 8.16 #From Azotobacter from #25 from 20 organized growth rates.xlsx#

RhRd = arange(0.1,10+0.1,0.1)   #(um/um) Radius of haptophyto / radius of diazotroph

Sf = 1-1/((RhRd)**3)#Space factor: see "Space taking factor computation.jpg"

MuC = Apho*(RhRd)**(3*B)/(C*(RhRd**(3*B)+1))
MuCs = Apho*(RhRd)**(3*B)*Sf/(C*(RhRd**(3*B)*Sf+1))
MuCs[RhRd<1] = nan

MuC2 = Apho2*(RhRd)**(3*B)/(C*(RhRd**(3*B)+1))
MuCs2 = Apho2*(RhRd)**(3*B)*Sf/(C*(RhRd**(3*B)*Sf+1))
MuCs2[RhRd<1] = nan

MuC3 = Apho3*(RhRd)**(3*B)/(C*(RhRd**(3*B)+1))
MuCs3 = Apho3*(RhRd)**(3*B)*Sf/(C*(RhRd**(3*B)*Sf+1))
MuCs3[RhRd<1] = nan

MuN = An2fix/(1+RhRd)**(3*B)
MuN2 = An2fix2/(1+RhRd)**(3*B)
MuN3 = An2fix3/(1+RhRd)**(3*B)
MuN4 = An2fix4/(1+RhRd)**(3*B)

Dpi = 300

figure(1)
plot(RhRd,MuC3,label='C lim Symbiosis')
plot(RhRd,MuN3,label='N lim Symbiosis', color="#D62728")
plot(RhRd,MuN2,label='N lim $\mathit{Crocosphaera}$')
plot(RhRd,MuN,label='N lim Heterocyst')

#plot(RhRd,MuCs,'--',label='C lim space')

plot(RhRd,MuN4,label='N lim Azotobacter',color='gray')
#plot(RhRd,MuC2,label='C lim (low)')
fill_between(RhRd, MuC, MuC2, where=MuC >= MuC2, facecolor='#DCE6F1',edgecolor = "none")

ylim(0,0.6)
xlim(left=0,right=10)
xlabel('R$_{H}$/R$_{D}$')
ylabel('$\mathit{\mu}$ (d$^{-1}$)')
legend(fontsize=17,loc=4,bbox_to_anchor=(1,0.2))
sf('','2 growth rates_multi',Dpi)

show()
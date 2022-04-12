'''
Created on Feb 17, 2020

@author: keiin
'''



from pylab import *
from E_BiosynthesisFromNH4 import *
from N2fixBudget import *
from FigSetting2 import *
from sf import *

rcParams['figure.figsize']=3.5,2.844
rcParams.update({'font.size': 9,
                 'lines.markersize':8,
                 'lines.markeredgewidth':1})
rcParams.update({'ytick.major.pad': 1})
rcParams.update({'xtick.major.pad': 1})
rcParams.update({'axes.linewidth':1}) 
rcParams.update({'lines.linewidth':1.5}) 

A = 0.216    #(various dimensions bassed on B) Factor converting V to Qc
B = 0.939    #(dimensionless) Power factor (Menden-Deuer and Lessard 2000)
Ycn = 6.3    #(molC molN-1) C:N
ep = 0.6
E = E_BiosynthesisFromNH4(ep) #(mol C mol C-1) CO2 production rate
Apho = (0.126+0.042)*12   #(d-1) Photosythesis rate divided by cellular C
Apho2 = (0.126-0.042)*12
Apho3 = 0.126*12

YcnN2fix,YresN2fix = N2fixBudget(ep) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
YresN2fix = 0
Y = (YcnN2fix+YresN2fix)/Ycn   #(energy cost) 
print(E,YresN2fix/Ycn)
print(YresN2fix)
C = (1+E+Y)       # Case without energy saving by light harvesting by UCYN-A
#C = (1+E+Y/2*1.5) # Case with energy saving by light harvesting by UCYN-A

#An2fix = 6.1/60*12 #(d-1) N2 fixation rate per N (based on 6.1 fmol h-1 cell-1 and fmol h-1 cell-1)
An2fix = 5*(18.845*4+61.0125)/61.0125   #value based on DDA
An2fix = 5.5*0.87*2.3272749437058176   #value based on DDA *5.5 times their own. 2.32... is the (Richelia N/heterocyst N) Value from D002 03 00

An2fix2 = 0.8  #for Crocosphaera 
An2fix3 = 10.42
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

MuN = An2fix/(1+RhRd**(3*B))
MuN2 = An2fix2/(1+RhRd**(3*B))
MuN3 = An2fix3/(1+RhRd**(3*B))
MuN4 = An2fix4/(1+RhRd**(3*B))

Dpi = 450

figure(2)
plot(RhRd,MuC3,label='C limited')
plot(RhRd,MuN3,label='N limited', color="#D62728")
#plot(RhRd,MuN2,label='N lim $\mathit{Crocosphaera}$')
#plot(RhRd,MuN,label='N lim DDA', color="#2CA02C")

#plot(RhRd,MuCs,'--',label='C lim space')

#plot(RhRd,MuN4,label='N lim Azotobacter',color='gray')
#plot(RhRd,MuC2,label='C lim (low)')
fill_between(RhRd, MuC, MuC2, where=MuC >= MuC2, facecolor='#DCE6F1',edgecolor = "none")

ylim(0,1)
xlim(left=0,right=10)
xlabel('R$_{H}$/R$_{D}$')
ylabel('$\mu$ (d$^{-1}$)')
ylim(0,1.38)
legend(loc=4,bbox_to_anchor=(1,0.2),edgecolor='#404040')
sf('','UCYN-A_only_Maranon',Dpi)

U = arange(size(RhRd))
for i in U:
    if MuC2[i] > MuN3[i]:
        RhRd_High = (RhRd[i] + RhRd[i-1])/2
        break

for i in U:
    if MuC3[i] > MuN3[i]:
        RhRd_Mid = (RhRd[i] + RhRd[i-1])/2
        break

for i in U:
    if MuC[i] > MuN3[i]:
        RhRd_Low = (RhRd[i] + RhRd[i-1])/2
        break
        
print(RhRd_Low,RhRd_Mid,RhRd_High)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# With Taniguchi data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
A = 0.216    #(various dimensions bassed on B) Factor converting V to Qc
B = 0.939    #(dimensionless) Power factor (Menden-Deuer and Lessard 2000)
Ycn = 6.3    #(molC molN-1) C:N
ep = 0.6
E = E_BiosynthesisFromNH4(ep) #(mol C mol C-1) CO2 production rate
Apho = 1.26*(1+E)     #(d-1) Photosythesis rate divided by cellular C
Apho2 = 0.76*(1+E)
Apho3 = 1.01*(1+E)

YcnN2fix,YresN2fix = N2fixBudget(ep) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
YresN2fix = 0
Y = (YcnN2fix+YresN2fix)/Ycn   #(energy cost) 
print(E,YresN2fix/Ycn)
print(YresN2fix)
C = (1+E+Y)       # Case without energy saving by light harvesting by UCYN-A
#C = (1+E+Y/2*1.5) # Case with energy saving by light harvesting by UCYN-A

#An2fix = 6.1/60*12 #(d-1) N2 fixation rate per N (based on 6.1 fmol h-1 cell-1 and fmol h-1 cell-1)
An2fix = 5*(18.845*4+61.0125)/61.0125   #value based on DDA
An2fix = 5.5*0.87*2.3272749437058176   #value based on DDA *5.5 times their own. 2.32... is the (Richelia N/heterocyst N) Value from D002 03 00

An2fix2 = 0.8  #for Crocosphaera 
An2fix3 = 10.42
An2fix4 = 8.16 #From Azotobacter from #25 from 20 organized growth rates.xlsx#

RhRd = arange(0.1,10+0.001,0.001)   #(um/um) Radius of haptophyto / radius of diazotroph

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

MuN = An2fix/(1+RhRd**(3*B))
MuN2 = An2fix2/(1+RhRd**(3*B))
MuN3 = An2fix3/(1+RhRd**(3*B))
MuN4 = An2fix4/(1+RhRd**(3*B))

figure(1)
plot(RhRd,MuC3,label='C limited')
plot(RhRd,MuN3,label='N limited', color="#D62728")
#plot(RhRd,MuN2,label='N lim $\mathit{Crocosphaera}$')
#plot(RhRd,MuN,label='N lim DDA', color="#2CA02C")

#plot(RhRd,MuCs,'--',label='C lim space')

#plot(RhRd,MuN4,label='N lim Azotobacter',color='gray')
#plot(RhRd,MuC2,label='C lim (low)')
fill_between(RhRd, MuC, MuC2, where=MuC >= MuC2, facecolor='#DCE6F1',edgecolor = "none")

ylim(0,1)
xlim(left=0,right=10)
xlabel('R$_{H}$/R$_{D}$')
ylabel('$\mu$ (d$^{-1}$)')
ylim(0,1.38)
legend(loc=4,bbox_to_anchor=(1,0.2),edgecolor='#404040')
sf('','UCYN-A_only_Taniguchi',Dpi)

U = arange(size(RhRd))

for i in U:
    if MuC2[i] > MuN3[i]:
        RhRd_High = (RhRd[i] + RhRd[i-1])/2
        break

for i in U:
    if MuC3[i] > MuN3[i]:
        RhRd_Mid = (RhRd[i] + RhRd[i-1])/2
        break

for i in U:
    if MuC[i] > MuN3[i]:
        RhRd_Low = (RhRd[i] + RhRd[i-1])/2
        break
        
print(RhRd_Low,RhRd_Mid,RhRd_High)

show()
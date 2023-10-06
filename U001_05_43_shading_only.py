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
                 'lines.markersize':6,
                 'lines.markeredgewidth':1})
rcParams.update({'ytick.major.pad': 1})
rcParams.update({'xtick.major.pad': 1})
rcParams.update({'axes.linewidth':1}) 
rcParams.update({'lines.linewidth':1.5}) 

A = 0.216    #(various dimensions bassed on B) Factor converting V to Qc
B = 0.866    #(dimensionless) Power factor (Non diatoms, Strathmann et al, 1967)
B = 0.863    #(dimensionless) Power factor (Verity et al, 1992)
B = 0.939    #(dimensionless) Power factor (Menden-Deuer and Lessard 2000)
Ycn = 6.3    #(molC molN-1) C:N
ep = 0.6
E = E_BiosynthesisFromNH4(ep) #(mol C mol C-1) CO2 production rate
Apho = 1.26*(1+E)     #(d-1) Photosythesis rate divided by cellular C
Apho2 = 0.76*(1+E)
Apho3 = 1.01*(1+E)
print(Apho,Apho2,Apho3)

YcnN2fix,YresN2fix = N2fixBudget(ep) #(mol C mol N-1) Converting N2 fixation to respiration and C consumption for electron
#YresN2fix = 0
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

Dpi = 450

figure(1)
plot(RhRd,MuC3,label='C limited default', color="#1f77b4")
plot(RhRd,MuC3+0.05,'--',label='C limited + 0.05 (d$^{-1}$)', color="#1f77b4")
plot(RhRd,MuC3+0.15,'-.',label='C limited + 0.15 (d$^{-1}$)', color="#1f77b4")
plot(RhRd,MuC3+0.3,':',label='C limited + 0.3 (d$^{-1}$)', color="#1f77b4")
plot(RhRd,MuN3,label='N limited default', color="#D62728")
plot(RhRd,MuN3+0.05,'--',label='N limited + 0.05 (d$^{-1}$)', color="#D62728")
plot(RhRd,MuN3+0.15,'-.',label='N limited + 0.15 (d$^{-1}$)', color="#D62728")
plot(RhRd,MuN3+0.3,':',label='N limited + 0.3 (d$^{-1}$)', color="#D62728")

ylim(0,1.2)
xlim(left=0,right=10)
xlabel('$\mathit{R_{H}/R_{D}}$')
ylabel('$\mathit{\mu}$ (d$^{-1}$)')


U = arange(size(RhRd))

def cp(MuC,MuN):  #Crossing point
    for i in U:
        if MuC[i] > MuN[i]:
            ic = i  #i crossing point
            break
    return ic

ic = cp(MuC3,MuN3)

ic05 = cp(MuC3+0.05,MuN3+0.05)
ic15 = cp(MuC3+0.15,MuN3+0.15)
ic3 = cp(MuC3+0.3,MuN3+0.3) #This sets one edge

ic0005 = cp(MuC3,MuN3+0.05)
ic0015 = cp(MuC3,MuN3+0.15)
ic003 = cp(MuC3,MuN3+0.30)

ic0500 = cp(MuC3+0.05,MuN3)
ic0515 = cp(MuC3+0.05,MuN3+0.15)
ic053 = cp(MuC3+0.05,MuN3+0.30)

ic1500 = cp(MuC3+0.15,MuN3)
ic1505 = cp(MuC3+0.15,MuN3+0.05)
ic153 = cp(MuC3+0.15,MuN3+0.30)

ic300 = cp(MuC3+0.3,MuN3) #This sets the other edge
ic305 = cp(MuC3+0.3,MuN3+0.05)
ic315 = cp(MuC3+0.3,MuN3+0.15)

# plot(RhRd[ic0005],MuC3[ic0005],'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic0015],MuC3[ic0015],'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic003],MuC3[ic003],'o',color='gray',markeredgecolor='gray')
#
# plot(RhRd[ic0500],MuC3[ic0500]+0.05,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic0515],MuC3[ic0515]+0.05,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic053],MuC3[ic053]+0.05,'o',color='gray',markeredgecolor='gray')
#
# plot(RhRd[ic1500],MuC3[ic1500]+0.15,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic1505],MuC3[ic1505]+0.15,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic153],MuC3[ic153]+0.15,'o',color='gray',markeredgecolor='gray')
#
# plot(RhRd[ic300],MuC3[ic300]+0.3,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic305],MuC3[ic305]+0.3,'o',color='gray',markeredgecolor='gray')
# plot(RhRd[ic315],MuC3[ic315]+0.3,'o',color='gray',markeredgecolor='gray')
#
# plot(RhRd[ic],MuC3[ic],'o',color='k')
# plot(RhRd[ic05],MuC3[ic05]+0.05,'o',color='k')
# plot(RhRd[ic15],MuC3[ic15]+0.15,'o',color='k')
# plot(RhRd[ic3],MuC3[ic3]+0.3,'o',color='k')

a = (RhRd[ic300],RhRd[ic003])
b = array([0,0])
c = array([2,2])
print(size(c >= b))

fill_between(a,b,c, where=c >= b, facecolor='#e6e6e6',edgecolor = "none")

#legend(loc=4,bbox_to_anchor=(1,0.2),edgecolor='#404040')
sf('','UCYN-A_sensitivity',Dpi)
#plot(RhRd[ic],MuN3[ic],'o', color='k')

figure(2)
plot([],[],label='C limited default', color="#1f77b4")
plot([],[],'--',label='C limited + 0.05 (d$^{-1}$)', color="#1f77b4")
plot([],[],'-.',label='C limited + 0.15 (d$^{-1}$)', color="#1f77b4")
plot([],[],':',label='C limited + 0.3 (d$^{-1}$)', color="#1f77b4")
plot([],[],label='N limited default', color="#D62728")
plot([],[],'--',label='N limited + 0.05 (d$^{-1}$)', color="#D62728")
plot([],[],'-.',label='N limited + 0.15 (d$^{-1}$)', color="#D62728")
plot([],[],':',label='N limited + 0.3 (d$^{-1}$)', color="#D62728")




legend(loc=4,bbox_to_anchor=(1,0.2),edgecolor='#404040')
sf('','UCYN-A_sensitivity_legend',Dpi)

show()
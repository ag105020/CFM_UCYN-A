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

VcVm = arange(0,10.01,0.1)
One = ones(size(VcVm))
PhRes = 1/6*(1+E)/E

MUcMUm = VcVm*PhRes*E/(1+E)

Dpi = 300

figure(2)
plot(VcVm,MUcMUm)
plot(VcVm,One,':',color='k')
xlabel('$\mathit{V_C/V_M}$')
ylabel('$\mathit{\mu_C/\mu_M}$')
xlim(0,10)
ylim(0,1.8)
sf('','MuCMuM',Dpi)

show()
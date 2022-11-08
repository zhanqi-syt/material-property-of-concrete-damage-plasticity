from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np

# concrete material parameters
ModelName = 'Model-1'                  # Add your own model name
MaterName = 'Concr-CDP'                # Add your own material name
fcr = 50.0                             # Compression strength
ecr = 0.000001*(700+172.*fcr**0.5)     # Peak strain under compression which is adjustable (empirical relationship)
e0 = 12285.0*fcr**0.264                # Young's modulus which is adjustable (empirical relationship)
miu = 0.2                              # Poisson's ratio
Sigma_y = fcr*0.4                      # Yield stress (0.4~0.6*fcr) which is adjustable (empirical relationship)
dens = 2.2e-09                         # Density
delte = 0.0001                         # Strain interval between input points
maxec = 0.01                           # Max compression strain
maxet = 0.002                          # Max tension strain

# Computing CDP curves
alpc = 0.157*fcr**0.785-0.905          # Declining parameter of stress-strain curves under compression (empirical relationship)
ftr = 0.1*fcr                          # Representative tension strength (empirical relationship)
etr = ftr**0.54*0.000065               # Stain for representative tension strength (empirical relationship)
alpt = 0.312*ftr**2                    # Declining parameter of stress-strain curves under tension (empirical relationship)
# Initialization
Epsilon_c = np.linspace(0.0, maxec, int(maxec/delte)+1)
Sigma_c = np.zeros(int(maxec/delte)+1)
Dam_c = np.zeros(int(maxec/delte)+1)
Epsilon_in_c = np.zeros(int(maxec/delte)+1)
Epsilon_pl_c = np.zeros(int(maxec/delte)+1)
n = e0 * ecr / (e0 * ecr - fcr)
rouc = fcr / (e0 * ecr)
for i in range(int(maxec/delte)+1):
    if Epsilon_c[i] <= ecr:
        Sigma_c[i] = Epsilon_c[i]*(rouc*n*e0/(n-1+(Epsilon_c[i]/ecr)**n))
    else:
        Sigma_c[i] = Epsilon_c[i]*(rouc*e0/(alpc*(Epsilon_c[i]/ecr-1)**2+Epsilon_c[i]/ecr))
    if Epsilon_c[i] != 0:
        Dam_c[i] = 1.0-np.sqrt(Sigma_c[i]/(Epsilon_c[i]*e0))
    Epsilon_in_c[i] = Epsilon_c[i]-Sigma_c[i]/e0
    Epsilon_pl_c[i] = Epsilon_c[i]-Sigma_c[i]/(e0*(1-Dam_c[i]))


Epsilon_t = np.linspace(0.0, maxet, int(maxet/delte)+1)
Sigma_t = np.zeros(int(maxet/delte)+1)
Dam_t = np.zeros(int(maxet/delte)+1)
Epsilon_in_t = np.zeros(int(maxet/delte)+1)
Epsilon_pl_t = np.zeros(int(maxet/delte)+1)
rout = ftr/(e0*etr)
for i in range(int(maxet/delte)+1):
    if Epsilon_t[i] <= etr:
        Sigma_t[i] = Epsilon_t[i]*e0*(rout*(1.2-0.2*((Epsilon_t[i]/etr)**5)))
    else:
        Sigma_t[i] = Epsilon_t[i]*e0*(rout/(alpt*((Epsilon_t[i]/etr-1.0)**1.7)+Epsilon_t[i]/etr))
    if Epsilon_t[i] != 0:
        Dam_t[i] = 1.0-np.sqrt(Sigma_t[i]/(Epsilon_t[i]*e0))
    Epsilon_in_t[i] = Epsilon_t[i]-Sigma_t[i]/e0
    Epsilon_pl_t[i] = Epsilon_t[i]-Sigma_t[i]/(e0*(1-Dam_t[i]))


start_c = np.argmax(Sigma_c > Sigma_y)
end_c = np.argmax(Dam_c > 0.9)
start_t = np.argmax(Sigma_t)
end_t = np.argmax(Dam_t > 0.9)
if end_c != 0:
    Harden_c = Epsilon_in_c[start_c:end_c].copy()
    Damage_c = Dam_c[start_c:end_c].copy()
    Harden_c[0] = 0.0
    Damage_c[0] = 0.0
    Harden_c = np.vstack((Sigma_c[start_c:end_c], Harden_c))
    Damage_c = np.vstack((Damage_c, Harden_c[1, :]))
else:
    Harden_c = Epsilon_in_c[start_c::].copy()
    Damage_c = Dam_c[start_c::].copy()
    Harden_c[0] = 0.0
    Damage_c[0] = 0.0
    Harden_c = np.vstack((Sigma_c[start_c::], Harden_c))
    Damage_c = np.vstack((Damage_c, Harden_c[1, :]))


if end_t != 0:
    Stiffen_t = Epsilon_in_t[start_t:end_t].copy()
    Damage_t = Dam_t[start_t:end_t].copy()
    Stiffen_t[0] = 0.0
    Damage_t[0] = 0.0
    Stiffen_t = np.vstack((Sigma_t[start_t:end_t], Stiffen_t))
    Damage_t = np.vstack((Damage_t, Stiffen_t[1, :]))
else:
    Stiffen_t = Epsilon_in_t[start_t::].copy()
    Damage_t = Dam_t[start_c::].copy()
    Stiffen_t[0] = 0.0
    Damage_t[0] = 0.0
    Stiffen_t = np.vstack((Sigma_t[start_t::], Stiffen_t))
    Damage_t = np.vstack((Damage_t, Stiffen_t[1, :]))


Harden_c = tuple([tuple(e) for e in Harden_c.T])
Damage_c = tuple([tuple(e) for e in Damage_c.T])
Stiffen_t = tuple([tuple(e) for e in Stiffen_t.T])
Damage_t = tuple([tuple(e) for e in Damage_t.T])

# Material Definitions
m = mdb.models[ModelName]
m.Material(name=MaterName)
mat = m.materials[MaterName]
mat.Density(table=((dens, ), ))
mat.Elastic(table=((e0, miu), ))
mat.ConcreteDamagedPlasticity(table=((30.0,0.1, 1.16, 0.667, 0.0001), ))
mat.concreteDamagedPlasticity.ConcreteCompressionHardening(table=Harden_c)
mat.concreteDamagedPlasticity.ConcreteTensionStiffening(table=Stiffen_t)
mat.concreteDamagedPlasticity.ConcreteCompressionDamage(table=Damage_c)
mat.concreteDamagedPlasticity.ConcreteTensionDamage(table=Damage_t)


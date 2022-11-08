from material import *
import numpy as np

# Concrete material parameters
ModelName = 'Model-1'                  # Add your own model name
MaterName = 'Concre'                   # Add your own material name
miu = 0.2                              # Poisson's ratio (from experimental data)
dens = 2.2e-09                         # Density (from experimental data)
fcr = 50.0                             # Compression strength (from experimental data)
e0 = 12285.0*fcr**0.264                # Elastic modulus (empirical relationship)(adjustable)
sy = fcr*0.4                           # Yield stress (0.4~0.6*fcr)(empirical relationship)(adjustable)
ecr = 0.000001*(700+172.*fcr**0.5)     # Strain for compression strength (empirical relationship)(adjustable)
ftr = 0.1*fcr                          # Tension strength (empirical relationship)(adjustable)
etr = ftr**0.54*0.000065               # Stain for tension strength (empirical relationship)(adjustable)
alpc = 0.157*fcr**0.785-0.905          # Descending params of compression curve (empirical relationship)(adjustable)
alpt = 0.312*ftr**2                    # Descending params of tension curve (empirical relationship)(adjustable)
dec = 0.0001                           # Strain increments in compression curve
dam_c_max = 0.9                        # Maximun compression damage
det = 0.0001                           # Strain increments in tension curve
dam_t_max = 0.9                        # Maximun tension damage

# Compression curve
n = e0 * ecr / (e0 * ecr - fcr)
rouc = fcr / (e0 * ecr)
cpr_all = []
dam_c = 0.0
sc = sy
ec = sy/e0
ec_in = 0
ec_pl = 0
i = 0
while dam_c < dam_c_max:
    i = i+1
    # cpr_all: compression information [No., strain, inelastic strain, plastic strain, stress, damage]
    cpr_all.append([i, ec, ec_in, ec_pl, sc, dam_c])
    ec = ec+dec
    if ec <= ecr:
        sc = ec*(rouc*n*e0/(n-1+(ec/ecr)**n))
    else:
        sc = ec*(rouc*e0/(alpc*(ec/ecr-1)**2+ec/ecr))
    dam_c = 1.0-np.sqrt(sc/(ec*e0))
    ec_in = ec - sc/e0
    ec_pl = ec - sc/(e0*(1-dam_c))

# Tension curve
rout = ftr/(e0*etr)
ten_all = []
dam_t = 0.0
st = ftr
et = etr
et_in = 0
et_pl = 0
i = 0
while dam_t < dam_t_max:
    i = i+1
    # ten_all: tension information [No., strain, inelastic strain, plastic strain, stress, damage]
    ten_all.append([i, et, et_in, et_pl, st, dam_t])
    et = et + det
    st = et*e0*(rout/(alpt*((et/etr-1.0)**1.7)+et/etr))
    dam_t = 1.0-np.sqrt(st/(et*e0))
    et_in = et-st/e0
    et_pl = et-st/(e0*(1-dam_t))

# Store in tuples
cpr_all = np.array(cpr_all)
ten_all = np.array(ten_all)
Harden_c = tuple([tuple(e) for e in cpr_all[:, [4, 2]]])
Damage_c = tuple([tuple(e) for e in cpr_all[:, [5, 2]]])
Stiffen_t = tuple([tuple(e) for e in ten_all[:, [4, 2]]])
Damage_t = tuple([tuple(e) for e in ten_all[:, [5, 2]]])

# Material Definitions
m = mdb.models[ModelName]
m.Material(name=MaterName)
mat = m.materials[MaterName]
mat.Density(table=((dens, ), ))
mat.Elastic(table=((e0, miu), ))
mat.ConcreteDamagedPlasticity(table=((30.0, 0.1, 1.16, 0.667, 0.0001), ))
mat.concreteDamagedPlasticity.ConcreteCompressionHardening(table=Harden_c)
mat.concreteDamagedPlasticity.ConcreteTensionStiffening(table=Stiffen_t)
mat.concreteDamagedPlasticity.ConcreteCompressionDamage(table=Damage_c)
mat.concreteDamagedPlasticity.ConcreteTensionDamage(table=Damage_t)

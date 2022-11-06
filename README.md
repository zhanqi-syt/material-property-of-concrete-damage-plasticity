If useful for your research, please cite: <br>
Yunhan Niu, Wenwu Wang, Yutai Su, Fengrui Jia, Xu Long. "Plastic damage prediction of concrete under compression based on deep learning." XXXXXX (2023): XXXXXX.<br>
Be free to use this code. If any questions, pls email us (suyutai@nwpu.edu.cn or xulong@nwpu.edu.cn). <br>

# material property of concrete damaged plasticity
## Keywords:
ABAQUS; Python script; CDP model; concrete damaged plasticity
## Introduction
Concrete damage plasticity (CDP) model is widely used to describe the elastoplastic properties of concrete. <br>
The tuples of CDP parameters in material definitions could be calculated by this code,
which could be used in Abaqus python script to automativally generate the material property of CDP. <br>
Demo is as follow: <br> <br> <br>
#Material Definitions <br>
m = mdb.models['Model-1'] <br>
m.Material(name='Concre') <br>
m.materials['Concre'].Density(table=((dens, ), )) <br>
m.materials['Concre'].Elastic(table=((e0, miu), )) <br>
m.materials['Concre'].ConcreteDamagedPlasticity(table=((30.0,0.1, 1.16, 0.667, 0.0001), )) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteCompressionHardening(table=Harden_c) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteTensionStiffening(table=Stiffen_t) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteCompressionDamage(table=Damage_c) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteTensionDamage(table=Damage_t) <br>

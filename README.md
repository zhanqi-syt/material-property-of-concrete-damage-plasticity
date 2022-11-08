If useful to your research, we would appreciate a citation:<br>
***Yunhan Niu, Wenwu Wang, Yutai Su, Fengrui Jia, Xu Long. 
"Plastic damage prediction of concrete under compression based on deep learning." 
XXXXXX (2023): XXXXXX.***<br>
Feel free to utilize this code. If any questions, please email us (suyutai@nwpu.edu.cn or xulong@nwpu.edu.cn). <br>

# Material property of concrete damaged plasticity in Abaqus python script
## Keywords:
ABAQUS; Python script; CDP model; concrete damaged plasticity
## Introduction
Concrete damage plasticity (CDP) model is widely used to describe the elastoplastic properties of concrete. 
The tuples of CDP parameters in material definitions could be calculated by this code,
which could be used in Abaqus python script to automativally generate the material property of CDP. <br>
Compression strength, density, and poisson's ratio are 3 necessary params, other params could be calculated according to empirical relationships and standards. <br>
Demo is as follow: <br> <br>
`# Material Definitions in Abaqus python script` <br>
m = mdb.models['Model-1'] <br>
m.Material(name='Concre') <br>
m.materials['Concre'].Density(table=((dens, ), )) <br>
m.materials['Concre'].Elastic(table=((e0, miu), )) <br>
m.materials['Concre'].ConcreteDamagedPlasticity(table=((30.0,0.1, 1.16, 0.667, 0.0001), )) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteCompressionHardening(table=Harden_c) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteTensionStiffening(table=Stiffen_t) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteCompressionDamage(table=Damage_c) <br>
m.materials['Concre'].concreteDamagedPlasticity.ConcreteTensionDamage(table=Damage_t) <br>

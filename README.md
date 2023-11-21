If useful to your research, we would appreciate a citation:<br>
***Yunhan Niu, Wenwu Wang, Yutai Su, Fengrui Jia, Xu Long. 
"Plastic damage prediction of concrete under compression based on deep learning." 
Acta Mechanica (Published: 20 October 2023)***<br>
Feel free to utilize this code. If any questions, please email us (suyutai@nwpu.edu.cn or xulong@nwpu.edu.cn). <br>

# Material property of concrete damage plasticity in Abaqus python script
## Keywords:
ABAQUS; Python script; CDP model; concrete damage plasticity
## Introduction
**Concrete damage plasticity (CDP) model** is widely used to describe the elastoplastic properties of concrete. <br><br>
***Compression strength***, ***density***, and ***Poisson's ratio*** are **ONLY 3 REQUIRED PARAMS**, 
and other params could be calculated according to empirical relationships and standards.<br><br>
The CDP tuples of material definitions are automatically generated by "main.py" in Abaqus python script. <br><br>
Abaqus interface is as follows: <br><br>
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

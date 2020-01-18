#######
# Modé́lisation d́e Réseaux Métaboliques
# MOCELL 2019
# Eden Darnige, Juan Manuel Garcia, Théo Gauvrit 
######

#Imports 
from cobra import *
import cobra.test
from cobra.flux_analysis import *

import matplotlib.pyplot as plt

##################################
#Une seconde analyse de flux (diapo 27)
##################################

##################################
#Analyse de la distribution des flux (diapo 28)
##################################

##################################
#Analyse de la variabilité́ des flux (diapo 29)
##################################

##################################
#Des KOs artificiels (diapo 30)
##################################

##################################
#Analyse de productibilité de molécules (diapos 31-32)
##################################
print("\nProductibilité de molécules\n")
S = Metabolite ( 
    id = 'S' , 
    formula='C6H10O12P2' , 
    name='S',
    compartment = 'c' )

T = Metabolite(
    id = 'T' , 
    formula='C3H5O6P' , 
    name='T',
    compartment = 'c' )

C1 = Metabolite(
    id = 'C1' , 
    formula='C3H5O6P' , 
    name='C1',
    compartment = 'c' )

C2 = Metabolite(
    id = 'C2' , 
    formula='C3H5O6P' , 
    name='C2',
    compartment = 'c' )

# S et T
rS = Reaction('rS')
rS.name = 'rS' 
rS.lower_bound = 0
rS.upper_bound = 10

rS.add_metabolites({ 
    S : 1.0
}) 

rT = Reaction('rT')
rT.name = 'rT' 
rT.lower_bound = 0
rT.upper_bound = 1000

rT.add_metabolites({ 
    T : -1.0
}) 

# Reactions 1-9 
r1 = Reaction('r1')
r1.name = 'r1' 
r1.lower_bound = 0
r1.upper_bound = 1000
r1.add_metabolites({ 
    S : -1.0,
    T : 1.0,
}) 

r2 = Reaction('r2')
r2.name = 'r2' 
r2.lower_bound = 0
r2.upper_bound = 1000
r2.add_metabolites({ 
    S : -1.0,
    C1 : -1.0,
    C2 : 1.0
}) 

r3 = Reaction('r3')
r3.name = 'r3' 
r3.lower_bound = 0
r3.upper_bound = 1000
r3.add_metabolites({ 
    C2 : -1.0,
    C1 : 1.0
}) 

r4 = Reaction('r4')
r4.name = 'r4' 
r4.lower_bound = 0
r4.upper_bound = 1000
r4.add_metabolites({ 
    C2 : -1.0,
    T : 1.0
}) 

r5 = Reaction('r5')
r5.name = 'r5' 
r5.lower_bound = 0
r5.upper_bound = 1000
r5.add_metabolites({ 
    S : -1.0,
    C1: -1.0,
    T : 1.0
}) 

r6 = Reaction('r6')
r6.name = 'r6' 
r6.lower_bound = 0
r6.upper_bound = 1000
r6.add_metabolites({ 
    C2 : -1.0,
    C1 : 2.0
}) 

r7 = Reaction('r7')
r7.name = 'r7' 
r7.lower_bound = 0
r7.upper_bound = 1000
r7.add_metabolites({ 
    T : 1.0,
    C1 : 1.0,
    C2 : -1.0
}) 

r8 = Reaction('r8')
r8.name = 'r8' 
r8.lower_bound = 0
r8.upper_bound = 1000
r8.add_metabolites({ 
    C1 : -1.0,
    S: -1.0,
    C2 : 2.0
}) 

r9 = Reaction('r9')
r9.name = 'r9' 
r9.lower_bound = 0
r9.upper_bound = 1000
r9.add_metabolites({ 
    C2 : -1.0,
    C1 : 1.0
}) 


#A
modelA = Model()

modelA = Model('ModelA')
modelA.add_reactions([rS,r1,rT])
modelA.objective="rT"

solution = modelA.optimize()
print("\nA - productible")
print(solution.fluxes)

#B
modelB = Model()

modelB = Model('modelB')
modelB.add_reactions([rS,r2,r3,r4,rT])
modelB.objective="rT"

solution = modelB.optimize()
print("\nB - non-productible")
print(solution.fluxes)

#C
modelC = Model()

modelC = Model('modelC')
modelC.add_reactions([rS,r5,rT])
modelC.objective="rT"

solution = modelC.optimize()
print("\nC - non-productible")
print(solution.fluxes)

#D
modelD = Model()

modelD = Model('modelD')
modelD.add_reactions([rS,r2,r4,r6,rT])
modelD.objective="rT"

solution = modelD.optimize()
print("\nD - productible")
print(solution.fluxes)

#E
modelE = Model()

modelE = Model('ModelE')
modelE.add_reactions([rS,r2,r7,rT])
modelE.objective="rT"

solution = modelE.optimize()
print("\nE - productible")
print(solution.fluxes)

#F
modelF = Model()

modelF = Model('modelF')
modelF.add_reactions([rS,r4,r8,r9,rT])
modelF.objective="rT"

solution = modelF.optimize()
print("\nF - productible")
print(solution.fluxes)

##################################
#Croissance sur différents substrats (diapo 35)
##################################
print("\nCroissance sur differents substrats")

#Import test model included in cobra
modelTextbook = cobra.test.create_test_model("textbook")

#Production of "ac" in function of O2 consumption 
prod_env = production_envelope(modelTextbook, 
    reactions = ["EX_o2_e"], objective="EX_ac_e",
    carbon_sources = "EX_glc__D_e")

prod_env2 = production_envelope(modelTextbook,
    ["EX_glc__D_e","EX_o2_e"])

print(prod_env)
print(prod_env2)

prod_env.plot(
    kind='line', x='EX_o2_e', y='flux_maximum');

plt.show()

##################################
#Des modèles plus gros (diapo 37)
##################################
modelEcoli = cobra.test.create_test_model("ecoli")
modelSalmonella = cobra.test.create_test_model("salmonella")


print("Le modèle a {} reactions".format(len(modelEcoli.reactions)))
print("Le modèle a {} reactions".format(len(modelSalmonella.reactions)))

print("Le modèle a {} génes".format(len(modelEcoli.genes)))
print("Le modèle a {} génes".format(len(modelSalmonella.genes)))

print("Le modèle a {} metabolites".format(len(modelEcoli.metabolites)))
print("Le modèle a {} metabolites".format(len(modelSalmonella.metabolites)))

# print (type(modelEcoli.reactions))

count=0
for i in modelEcoli.reactions: 
    if i in modelSalmonella.reactions: 
        count+=1
print ("There are {} reactions in common".format(count))

count=0
for i in modelEcoli.metabolites: 
    if i in modelSalmonella.metabolites: 
        count+=1
print ("There are {} metabolites in common".format(count))


def findBiomassMetabolites(model): 
    submodel=model.reactions.query("biomass",attribute='id')
    if (submodel==[]):
        submodel=model.reactions.query("BIOMASS",attribute='id')
    for i in submodel: 
        print(i)
    met=[]  
    for i in submodel: 
        for j in i.metabolites: 
            if j.name not in met: 
                met.append(j.name)
    return met 

# print("EColi: \n")
reac1=findBiomassMetabolites(modelEcoli)
# print(reac1)
print("Salmonella: \n")
reac2=findBiomassMetabolites(modelSalmonella)
# print(reac2)
count=0
for j in reac1: 
    if j in reac2: 
        count+=1
print("\nAu total, il y a {} metabolites en commun qui font partie de la biomasse".format(count)) 

## FBA Analysis 
# double_del_eColi=double_gene_deletion(modelEcoli)
# print(list(del_eColi.index.values))
# print(del_eColi.summary())
# print(del_eColi.index.values[del_eColi.index(max(del_eColi.growth))])
sol1_eco=modelEcoli.optimize()

sol2_salmo=modelSalmonella.optimize()
print(sol1_eco.objective_value)
print(sol2_salmo.objective_value)

def setReacto0(model,react_name): 
    medium=model.medium
    medium[react_name]=0.0
    model.medium=medium
    return model

## Anaerobic condition 
ana_Ecoli=setReacto0(modelEcoli,"EX_o2_e")
ana_Salmo=setReacto0(modelSalmonella,"EX_o2_e")

sol1_ana=ana_Ecoli.optimize()
sol2_ana=ana_Salmo.optimize()

print(sol1_ana.objective_value)
print(sol2_ana.objective_value)

## Carbon input limited
noC_Ecoli=setReacto0(modelEcoli,"EX_co2_e")
noC_Ecoli=setReacto0(noC_Ecoli,"EX_glc__D_e")

noC_Salmo=setReacto0(modelSalmonella,"EX_co2_e")
noC_Salmo=setReacto0(noC_Salmo,"EX_glc__D_e")

sol1_noC=noC_Ecoli.optimize()
sol2_noC=noC_Salmo.optimize()

print(sol1_noC.objective_value)
print(sol2_noC.objective_value) 
##################################
#Un peu de bioengineering (diapo 38)
##################################
modelEcoli=cobra.test.create_test_model("ecoli")

amiadi=Metabolite(id="amndp",formula="C6H11NO4",name="L_2_aminoadipate",compartment="c")
peni_int=Metabolite(id="am_carb_cys_val",name="N-[(5S)-5-amino-5-carboxypentanoyl]-L-cysteinyl-D-valine",formula="C14H24N3O6S",compartment="c")
isopen_N=Metabolite(id="iso_N",name="Isopenicilin_N",formula="C14H20N3O6S",compartment="c")

val=modelEcoli.metabolites.get_by_id("val__L_c")
cys=modelEcoli.metabolites.get_by_id("cys__L_c")

adp_entry=Reaction(
    id="adp_entry",
    name="Aminoadipate entry",
    upper_bound=1000,
    lower_bound=0
)

adp_cys_val=Reaction(
    id="adp_cys_val_synth",
    name="δ-(L-α-aminoadipyl)-L-cysteinyl-D-valine synthetase",
    lower_bound=0,
    upper_bound=1000
)
isopen_synthase=Reaction(
    id="isopen_synth",
    name="Isopenicilin-N synthase",
    upper_bound=1000,
    lower_bound=0
)
isopen_prod=Reaction(
    id="iso_prod",
    name="Isopenicilin production",
    upper_bound=1000,
    lower_bound=0
)
adp_entry.add_metabolites({
    amiadi: 1.0
})

adp_cys_val.add_metabolites({
    amiadi: -1.0,
    val: -1.0,
    cys: -1.0,
    modelEcoli.metabolites.h2o_c: -1.0,
    modelEcoli.metabolites.atp_c: -3.0,
    modelEcoli.metabolites.ppi_c: 3.0,
    modelEcoli.metabolites.amp_c: 3.0,
    modelEcoli.metabolites.h_c: 3.0,
    peni_int: 1.0

})
isopen_synthase.add_metabolites({
   peni_int: -1.0,
   modelEcoli.metabolites.o2_c: -1.0,
   modelEcoli.metabolites.h2o_c: 2.0,
   isopen_N: 1.0 
})
isopen_prod.add_metabolites({
    isopen_N: 1.0
})

EColi_trans=modelEcoli
EColi_trans.add_reactions([adp_entry,isopen_prod, adp_cys_val,isopen_synthase])

medium=EColi_trans.medium
medium["EX_glc__D_e"]=1000.0
EColi_trans.medium=medium

sol_trans=EColi_trans.optimize()

print(sol_trans.objective_value)
print(" isopeniciline production: " + str(sol_trans.fluxes["iso_prod"]))

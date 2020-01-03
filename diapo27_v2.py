## Modele simple 2
reac1=Reaction("Reaction1")
reac1.name="Reaction 1"
reac1.lower_bound=0
reac1.upper_bound=1000
reac1.add_metabolites({
    A: -1.0,
    B: -2.0,
    C: 4.0
})

reac2=Reaction("Reaction2")
reac2.name="Reaction 2"
reac2.upper_bound=1000
reac2.lower_bound=0
reac2.add_metabolites({
    A: -2.0,
    D: -1.0,
    E: -1.0,
    F: 2.0
})

reac3=Reaction("Reaction3")
reac3.name="Reaction 3"
reac3.upper_bound=1000
reac3.lower_bound=0
reac3.add_metabolites({
    F: -1.0,
    C: 1.0
})

reac4=Reaction("Reaction4")
reac4.name="Reaction 4"
reac4.upper_bound=1000
reac4.lower_bound=0
reac4.add_metabolites({
    B: -1.0,
    D: 1.0,
    E: 1.0
})
r_imp=Reaction("Input")
r_imp.name="Input reaction"
r_imp.upper_bound=10
r_imp.lower_bound=0
r_imp.add_metabolites({
    B: 1.0
})

r_bio=Reaction("Biomass")
r_bio.name="Biomass synthetysation"
r_bio.upper_bound=1000
r_bio.lower_bound=0
r_bio.add_metabolites({
    C: -2.0,
    Biomass: 1.0

})
r_im_A=Reaction("A_input")
r_im_A.name="Input of A"
r_im_A.upper_bound=20
r_im_A.lower_bound=0
r_im_A.add_metabolites({
    A: 1.0
})

model_simple2=Model("Ex2_modele")
model_simple2.add_reactions([reac1,reac2,reac3,reac4,r_imp,r_bio,r_exp_bio,r_im_A])

model_simple2.objective="Expression_bio"

for reac in model_simple2.boundary: 
    print (reac)
solution=model_simple2.optimize()
for reac in model_simple2.reactions: 
    print(str(reac.id) + " : " + str(solution.fluxes[reac.id]))

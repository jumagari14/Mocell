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


#Une seconde analyse de flux (diapo 27)

#Analyse de la distribution des flux (diapo 28)

#Analyse de la variabilité́ des flux (diapo 29)

#Des KOs artificiels (diapo 30)

#Analyse de productibilité de molécules (diapos 31-32)

#Croissance sur différents substrats (diapo 35)
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
    kind='line', x='EX_o2_e', y='carbon_yield_maximum');

plt.show()

#Des modèles plus gros (diapo 37)

#Un peu de bioengineering (diapo 38)

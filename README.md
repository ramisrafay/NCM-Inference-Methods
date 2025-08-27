# How to quantify immigration from community abundance data using the Neutral Community Model 
This is companion code for the manuscript "How to quantify immigration from community abundance data using the Neutral Community Model", by Ramis Rafay, Eric Jones, David Sivak, and Jane Fowler.

# Steady state simulations 
Code for fitting the neutral community models to abundance/sequence datasets and inferring immigration rates are available as R scripts. Running "Paper_simulations_figures.R" reproduces figures 1-5 in the main manuscript (Ensure that "Function_definitions.R" is also present in your working directory). 

R code written by Ramis Rafay (ramisrafay@gmail.com) -- email if you have
any questions.

# Dynamic simulations 
The file dynamic_NCM_simulations_cleaned_EJ.py contains python code (version 3.13.5) to perform dynamic simulations of the Neutral Community Model (NCM).

This program supports simulations based on the classic mainland/island NCM (i.e., with source/local communities) and for a metapopulation NCM in which the aggregate of all other islands acts as a mainland for a focal island.  This code reproduces Figures S1 and S12.

By default, the code runs a Minimal Working Example (MWE) of a dynamic NCM simulation. Modify lines at the end of file to modify this MWE. To generate Figures S1 and S12, uncomment the corresponding lines at the end of the file, but note that these figures take many hours to run by default. 

Python code written by Eric Jones (jones.eric93@gmail.com) -- email if you have
any questions.

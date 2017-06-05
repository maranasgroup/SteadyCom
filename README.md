# SteadyCom
Computing the maximum growth rate of a microbial community at community steady-state.

Chan SHJ, Simons MN, Maranas CD (2017) SteadyCom: Predicting microbial abundances while ensuring community stability. PLoS Comput Biol 13(5): e1005539. https://doi.org/10.1371/journal.pcbi.1005539

Main functions:

createCommModel: create a community model from individual models

SteadyComCplex: given a community model constructed by createCommModel, find the maximum community growth rate by the SteadyCom procedure

SteadyComFVACplex: perform FVA under the SteadyCom framework

SteadyComPOACplex: perform pairwise FVA under the SteadyCom framework

See SteadyCom/doc/index.html for detailed documentations

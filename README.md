# SteadyCom
Computing the maximum growth rate of a microbial community at community steady-state.
Ref:
Chan SHJ, Simons MN, Maranas CD. (Accepted) SteadyCom: Predicting Microbial Abundances while Ensuring Community Stability. PLOS Comput Biol.

Main functions:

createCommModel: create a community model from individual models

SteadyComCplex: given a community model constructed by createCommModel, find the maximum community growth rate by the SteadyCom procedure

SteadyComFVACplex: perform FVA under the SteadyCom framework

SteadyComPOACplex: perform pairwise FVA under the SteadyCom framework

See SteadyCom/doc/index.html for detailed documentations

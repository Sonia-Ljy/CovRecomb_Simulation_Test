# CovRecomb_Simulation_Test
The simulation tests for the CovRecomb method by simulating SARS-CoV-2 evolution.


## Overview
A forward-evolution simulator called [CovSimulator](https://github.com/weigangq/cov-db/blob/master/scripts/CovSimulator.py) (Saymon Akther, 2021) was used to generate the simulation recombination dataset by considering SARS-CoV-2 genome evolution and transmission regular. A few changes were made to CovSimulator to meet the requirements of our simulation test.

### Two stages in simulation
The simulation was divided into two stages: (1) the initial lineage generation process and (2) the lineage evolution process. In the first stage, starting from the SARS CoV-2 Wuhan strain, the genome would experience n generations (G) to generate 2^n sequences. Then, we randomly sampled several sequences in the last generation and took them as the initial composition of viral lineages. In the second stage, apart from the similar evolution process as stage one, homologous inter-lineage and intra-lineage recombination would occur with a Poisson distributed rate of population size for each generation. We ensured that at least one inter-lineage recombinant and one intra-lineage recombinant occurred at each generation regardless of the population size. Notably, during each generation in the two stages, all genomes would mutate at random positions with the number of sites following the Poisson distribution. A preset proportion of sequences would share a preset number of homologous mutations, which would help to integrate the convergent evolution factor into the simulation.

## Workflow
<img src="img/workflow.png"/>


## Acknowledgements
We would like to express our special thanks to the [CovSimulator](https://github.com/weigangq/cov-db/blob/master/scripts/CovSimulator.py) developers for their contributions in the field of simulation testing.

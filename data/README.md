## Benchmark Instances

Except for the **NewHard** instances and their instance generator (file instance_generator_NewHard.py in this directory) for the BPP and CSP proposed in this work, all other instances were taken from previous studies in the literature. All instances are available in this repository and follow the **BPP format**, meaning that items are not grouped by size and have unit demand. Below, we detail the source of each set of instances.

### BPP and CSP Instances

These instances are part of the [BPPLib](https://site.unibo.it/operations-research/en/research/bpplib-a-bin-packing-problem-library), a comprehensive library of challenging instances for the BPP and the CSP. For more information, refer to:

> Delorme, M., Iori, M., & Martello, S. (2018). **BPPLIB: A library for bin packing and cutting stock problems**. *Optimization Letters*, 12(2), 235–250.

### SSP Instances

The SSP instances are organized into three categories:

- CSP instances from BPPLib  
- A1, A2, A3, and B instances  
- GI instances  

The A and B instance classes are available at [https://github.com/wotzlaff/ssp-data](https://github.com/wotzlaff/ssp-data) and were originally introduced in:

> Martinovic, J., Delorme, M., Iori, M., Scheithauer, G., & Strasdat, N. (2020). **Improved flow-based formulations for the skiving stock problem**. *Computers & Operations Research*, 113, 104770.

**Note:** Although class A3 instances were generated in this work, they were only used in subsequent research.

The GI instances consist of CSP instances with varying numbers of item types (125, 250, 500, 750, and 1000). They were proposed in:

> Gschwind, T., & Irnich, S. (2016). **Dual Inequalities for Stabilized Column Generation Revisited**. *INFORMS Journal on Computing*, 28(1), 175–194.

However, only the 125, 250, and 500-instance classes were used in that original work and are partially available in BPPLib. The complete set—including 750 and 1000 item types—was later used in SSP studies and can be found at:  
[https://logistik.bwl.uni-mainz.de/forschung/benchmarks/](https://logistik.bwl.uni-mainz.de/forschung/benchmarks/)

### IPMS Instances

These instances were proposed in:

> Mrad, M., & Souayah, N. (2018). **An arc-flow model for the makespan minimization problem on identical parallel machines**. *IEEE Access*, 6, 5300–5307.

We obtained access to them via the authors of subsequent work:

> Gharbi, A., & Bamatraf, K. (2022). **An improved arc flow model with enhanced bounds for minimizing the makespan in identical parallel machine scheduling**. *Processes*, 10(11).

### CCBPP Instances

This instance set was kindly provided by the original authors upon request. They were introduced in:

> Borges, Y.G., Miyazawa, F.K., Schouery, R.C., & Xavier, E.C. (2020). **Exact algorithms for class-constrained packing problems**. *Computers & Industrial Engineering*, 144, 106455.

### OOEBPP Instances

The first part of the instances comes from [2DPackLIB](https://site.unibo.it/operations-research/en/research/2dpacklib), a benchmark library for two-dimensional orthogonal cutting and packing problems. The corresponding publication is:

> Iori, M., de Lima, V.L., Martello, S., & Monaci, M. (2021). **2DPackLib: A two-dimensional cutting and packing library**. *Optimization Letters*.

The second part consists of an adaptation of the **Random Class** BPP instances to the OOEBPP context, as proposed by:

> de Lima, V.L., Iori, M., & Miyazawa, F.K. (2023). **Exact solution of network flow models with strong relaxations**. *Mathematical Programming*, 197, 813–846.

This second set was kindly provided by the original authors upon request.
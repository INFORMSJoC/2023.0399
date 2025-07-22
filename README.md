[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Solving Cutting Stock Problems via an Extended Ryan-Foster Branching Scheme and Fast Column Generation

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [CC BY-NC 4.0 License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper
[Solving Cutting Stock Problems via an Extended Ryan-Foster Branching Scheme and Fast Column Generation](https://doi.org/10.1287/ijoc.2023.0399) by Renan F. F. da Silva and Rafael C. S. Schouery.


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0399

https://doi.org/10.1287/ijoc.2023.0399.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Silva2025cd,
  author    = {Silva, Renan F. F. da and Schouery, Rafael C. S.},
  publisher = {INFORMS Journal on Computing},
  title     = {Solving Cutting Stock Problems via an Extended Ryan-Foster Branching Scheme and Fast Column Generation},
  year      = {2025},
  doi       = {10.1287/ijoc.2023.0399.cd},
  url       = {https://github.com/INFORMSJoC/2023.0399},
  note      = {Available for download at https://github.com/INFORMSJoC/2023.0399}
}
```

---

## Description

In the paper, we present a general framework to address several cutting and packing problems, including:

- Cutting Stock Problem (CSP)
- Skiving Stock Problem (SSP)
- Identical Parallel Machines Scheduling with Minimum Makespan (IPMS)
- Ordered Open-End Bin Packing Problem (OOEBPP)
- Class-Constrained Bin Packing Problem (CCBPP)

---

## Requirements

In our experiments, all solvers were compiled using **GCC 11.4.0** and executed on **Ubuntu 22.04** with **Gurobi 10.0.1**.
By default, Gurobi is assumed to be installed in `/opt/gurobi1001/linux64`.

> If you are using a different installation path or a different version of Gurobi, please update the corresponding `Makefile` located at `src/<solver-name>/Makefile`.

---

## Running the Solver

You may run either a single instance or a complete set of instances. In our experiments, we used the Python scripts located in `scripts`, which:

1. Compile the solver and the solution checker
2. Execute the solver on each instance
3. Generate output files in `src/<solver-name>/out`:
   - `instance-name.gb`: Gurobi log file
   - `instance-name.log`: computed metrics
   - `instance-name.sol`: solution file

These scripts handle compilation and execution automatically for all instances in the set.

---

### Running a Single Instance Manually

To run a single instance manually, follow these steps:

```
cd src
make checkerCSP
cd CSPSolver
make
./CSPSolver ../../instances/CSP/AI202/201_2500_DI_0.txt ../checkerCSP
```

---

## Paper's Computational Experiments

The full set of computational results presented in the paper can be found in the `results/` directory.

## Acknowledgment

This research was supported by the São Paulo Research Foundation (FAPESP) grants #2015/11937-9 and #2022/05803-3;
the Brazilian National Council For Scientific and Technological Development (CNPq) grant #312345/2023-2;
and the Teaching, Research and Extension Support Fund of University of Campinas (FAEPEX/UNICAMP) grant #2372/23.

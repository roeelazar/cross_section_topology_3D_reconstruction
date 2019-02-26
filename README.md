# Robust optimization for topological surface reconstruction

<p align="center">
  <img src="silicon.png" width="448">
</p>

This repository contains the source code for the SIGGRAPH 2018 paper: [Robust optimization for topological surface reconstruction](https://dl.acm.org/citation.cfm?id=3201348).

This paper proposes an optimization method for surface reconstruction under topological constraints. The input to our method is a prescribed genus for the reconstructed surface, a partition of the ambient volume into cells, and a set of possible surface candidates and their associated energy within each cell. The outputs are one candidate per cell so that their union is a connected surface with the prescribed genus that minimizes the total energy.
Beside the core engine, an example is attached.

## Getting Started

### Setting up the environment
The code is compatible with MATLAB 2018a.

#### Additional tools
In order to run the code, one should download [Yalmip optimization package](https://yalmip.github.io/download/) and add it to the path.
Next, a solver that is compatible with yalmip and semidefinite programming should be installed. In the silicon example and in the paper mosek (free for academic use) was used (information about other solvers can be found [here](https://yalmip.github.io/allsolvers/)).

#### Mex files
For optimization running time reduction purposes, a mex file was generated. If you are using a Windows 64 bit, then no further action is needed. Otherwise you need to mex the cpp file located in Utils\mincutMex.cpp.

### Downloads
Clone the repository and work from the project directory. Add to the path the folders Utils, Core and Data_structure.

## Running 
### Silicon example
After all the above was executed, type in the command window:
```
silicon_example
```

### Preparing a reconstruction problem
To generate a new problem, first formulate your problem as mentioned in the documentaion of Core\BnBOptimization,
then generate BnBroot type and insert it to BnBOptimization(BnBroot).

## Running time note
The optimization's running time has been significantly decreased since the paper's publication.

## Acknowledgments
This paper was supported by the European Research Council, ERC-cog grant no. 771136-LiftMatch,
and the Israel Science Foundation, grant no. ISF 1830/17. ZH and
TJ acknowledge the support of NSF grants IIS-0846072, IIS-1302200,
RI-1618685.

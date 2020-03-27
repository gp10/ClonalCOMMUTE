# ClonalCOMMUTE :: Modeling Clonal COMpetition dynamics in MUTated Epithelium
Computational methods to explore cell competition and clone dynamics in mutated mouse epithelium.

This repository contains code used to simulate epithelial clone dynamics under the Neighbor-Constrained Fitness (NCF)-model hypothesis, as defined in the publication:
  > Colom B, Alcolea MP, Piedrafita G, Wabik A, Hall MWJ, Dentro SC, Fowler JC, Herms A, King C, Ong SH, Sood RK, Gerstung M, Martincorena I, Hall BA, Jones PH (2020) Spatial competition shapes the dynamic mutational landscape of normal esophageal epithelium. _accepted in Nat. Genet._

### Graphical abstract
![GraphicalAbstract](https://github.com/gp10/ClonalCOMMUTE/blob/master/Graphical_abstract_ClonalCOMMUTE.png)

### Overview of the NCF-model and implementation
The **NCF model** grounds on the plasticity of the stochastic progenitor cell fates within the epithelium, cells which debate between division and terminal differentiation, and assumes that mutant cells battle for a limited space with different resilience to differentiate, and fitness advantages that are not cell-intrinsic but relative to cell neighborhood. In the manuscript we show how this model accommodates experimental data and has major implications for clonal growth, competition and selection.

The NCF model is simulated by means of a spatial, lattice-based stochastic Monte Carlo implementation.

In this implementation, the basal proliferative compartment of epithelium is represented by a 2D grid where every position is occupied by a single progenitor cell.
- In the **Simple** `ProtocolType` scenario, there is initially a mixed population of normal WT cells (N) and Mutant cells (M) induced by treatment with mutagen DEN (with induction rate X<sub>mut</sub>) (t0). N and M cells show the same division potential, and once a cell is randomly chosen to divide it replaces a neighbor (replacement rate: parameter Lambda). However, M and N populations differ in their resilience to be displaced (differentiation) by other cells (parameter delta<sup>M</sup>). Individual cell fates are determined by their value of delta<sup>M</sup> relative to that of their neighbors. As a consequence mutant fitness gain, and hence clonal growth and selection, is not an intrinsic, cell-autonomous property but sensitive to spatial constraints.

Furthermore, two more complicated scenarios can be simulated and used to challenge our model inferences:
- The **DEN_IndMml** `ProtocolType` recreates a protocol where initial mutagenesis (by DEN; induction rate: X<sub>mut</sub>, mutant fitness: delta<sup>M</sup>) is followed by a late induction (rate: freqLabel) of highly competitive mutant cells (delta<sup>Mml</sup>) after a period of chase (LagTime).
- The **IndMml_Ind** `ProtocolType` recreates a protocol where there is an initial induction (rate: freqLabel) of highly competitive mutant cells (delta<sup>Mml</sup>) followed by induction (freqLabel2) of random (Confetti) labelling after a period of chase (LagTime).

### Model parameters and user options

#### Protocol choice
With the tools presented here the user can choose (`selectProtocol`) between the different protocol scenarios (`ProtocolType`) introduced above, each offering slightly different parameter controls. In turn, in the main **Simple** scenario, two versions of the model exist:
- **monoMut** : Here it is assumed for simplicity that all mutant cells M take the same value of fitness (delta<sup>M</sup>), i.e. they share the same mutated genotype.
- **polyMut** : Here every single mutant clone M<sub>1</sub>, M<sub>2</sub>, .., M<sub>m</sub> is assigned a different fitness value, with delta<sup>M</sup> randomly drawn from a distribution F=(1-Gamma(κ,1/κ)), with shape determined by parameter κ.
Mutagenesis in the more complicated scenarios is simulated as in **polyMut**.

#### List of adjustable NCF-model parameters
| Param.     | Description |
| --------   | ----------- |
| Lambda     | Fixed rate of cell replacement (same for N and M cells, consistent with _in vivo_ proliferation assays). |
| X<sub>mut</sub>     | Rate (%) of mutant cell induction caused by DEN at t=0. |
| delta<sup>M</sup>   | Mutant M fitness propensity: M resilience to differentiate (used in **monoMut** case only). |
| κ          | Shape parameter of the fitness landscape, where each mutant M~i~ fitness propensity delta<sup>Mi</sup> is randomly drawn from (used in **polyMut** case only).  |
| delta<sup>Mml</sup> | Mml mutant fitness propensity: Mml resilience to differentiate (used in **DEN_IndMml** and **IndMml_Ind** settings only). |

#### Other adjustable simulation parameters
| Param.        | Description |
| --------      | ----------- |
| freqLabel 	| Rate (%) of cells labelled by a fluorescent reporter (or Mml, in the case of **DEN_IndMml** and **IndMml_Ind**). |
| freqLabel2    | Rate (%) of cells labelled by a fluorescent reporter at late protocol stage (in the case of **IndMml_Ind** setting). |
| lattice.Dim   | 2D grid size (number of cells per dimension of the square lattice). |
| lattice.Neigh | Neighborhood geometry (4, 6 or 8 neighbors per cell). |
| timelim       | Simulation time span post-labelling (weeks) (DEN treatment is considered instantaneous at t=0). |
| nval        	| Number of (regularly spaced) time points when variables are retrieved and saved. |

### Main scripts
- **Analysis-Sim2DCompetition-BasalCloneSizes.m** : main script to run 2D simulations of clone competition and plot results.
  + The user can select the protocol (`selectProtocol`) and parameter values of interest and simulations are run (a simulation is included under same parameter conditions except that no mutagenesis occurs at time=0 as a CTL).
  + Results are saved in nested `Datasets` folder.
  + Plots are generated for the following clonal features:
    * average clone size over time
    * clone density over time
    * basal clone size distributions over time
    * 2D views of epithelial layer highlighting the clones over time
    * total fraction of basal mutant cells over time
    * fraction of (fluorescently)-labelled basal cells over time
- **Launch-Sim2DCompetition-examples.m** : script that launches example simulation cases under the different protocol scenarios, intended to help understanding the syntax and tuning parameter values under the protocol of interest.

### Dependencies
- MonteCarloSimulator-2Dgrid-SP-MutCloneDynamics.m : main function called to run lattice-based simulations of basal progenitor cell competition in CTL or DEN conditions.
- MonteCarloSimulator-2Dgrid-SP-MutCloneDynamics-challenge.m : main function called to run lattice-based simulations of basal progenitor cell competition in the more complex, two-stage 'challenge' protocols involving Mml* induction post-DEN or random labelling post-Mml* induction.
- SelectModelParamVal.m : function called to load specific, preset model parameter values.
- calculate-AvgCloneSize.m : function called to calculate the average clone size from clone size distributions at different time points.
- calculate-CloneDens.m : function called to calculate the total % of surviving clones (clone density) at different time points.
- calculate-CloneSizeDist.m : function called to calculate the normalized, cumulative clone size distributions at different time points.
- plot-AvgCloneSize.m : function called to plot the average clone size over time.
- plot-CloneDens.m : function called to plot the total % of surviving clones (clone density) over time.
- plot-CloneSizeDist.m : function called to plot the normalized, cumulative clone size distribution at given time points.
- plot-2Dview.m : function called to plot a spatial 2D view of clones (labelled in different colors) at given time points.
- plot-FreqMut.m : function called to plot the overall fraction (%) of mutant cells over time.
- plot-FreqLabel.m : function called to plot the fraction (%) of fluorescent labelled cells over time.
- subsampling-clones-overtime.m : function called to build random subsets of simulated clones (sampled by random permutation from the pool of simulated data) for the purpose of plausible intervals calculations.
- `Datasets` folder : contains some example, preset datasets used in the manuscript and stores new data generated by user.

### Requirements
Matlab R2016b


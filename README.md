# Visualize Binary Trees

### **Kyle Cranmer, Sebastian Macaluso and Duccio Pappadopulo**

Note that this is an early development version. 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 

## Introduction

This is a visualization package for binary trees. In particular, this packages allows to compare visualizations of jets generated with the [`Toy Generative Model for Jets`](https://github.com/SebastianMacaluso/ToyJetsShower) package. 

There are three main notebooks:

1) [`binaryTreeStructure`](binaryTreeStructure.ipynb): this  notebook explains the structure of the jet tree dictionaries (generated with the [`Toy Generative Model`](https://github.com/SebastianMacaluso/ToyJetsShower)  and some examples on how to access the information contained in the dictionaries. Also, any binary tree that satisfies the structure described in this notebook, could be visualized with this package.


2) [`visualize1DTree.ipynb`](visualize1DTree.ipynb): contains the *1D tree-only* jet visualizations.


3)[`visualizeClusterMap.ipynb`](visualizeClusterMap.ipynb): contains the *2D heat clustermaps* of the jets.


## 2D heat clustermaps


Given a tree structure, 2D heat clustermaps visualizations are obtained with [`visualizeClusterMap.ipynb`](visualizeClusterMap.ipynb).
Specifically, given a jet (the set of its constituents), visualization for the different trees obtained with each clustering algorithm (and the truth tree) are displayed. 

There are 2 functions for visualization:

a) function heat_dendrogram: Create  a heat dendrogram clustermap. If two jets are given as input:
- 1) If only one jet is given as input, both rows and columns are ordered according to that jet tree.
- 2) If truthJet and recluster_jet1, the visualization shows the heat data map of the truth jet, with columns ordered according to the truth jet and rows according to recluster_jet1.
- 3) If recluster_jet1 and recluster_jet2, the visualization shows the heat data map of recluster_jet2, with columns ordered according to the recluster_jet1 and rows according to recluster_jet2.

b) function dendrogramDiff:
Given two jet algorithms heat data matrices, we reorder the heat matrices according to the truth jet order (order in which leaves are accessed when traversing the truth tree) and take the difference.

Both functions have the same arguments:
- Truth jet dictionary
- recluster_jet1: reclustered jet 1
- recluster_jet2: reclustered jet 2
- full_path: Bool. If True, then use the total number of steps to connect a pair of leaves as the heat data. If False, then Given a pair of jet constituents {i,j} and the number of steps needed for each constituent to reach their closest common ancestor {Si,Sj}, the heat map scale represents the maximum number of steps, i.e. max{Si,Sj}.
- FigName: Dir and location to save a plot.


Usage:
- Input_dir: select the directory with the jets dictionaries.
- Input_jet: select the jet filename for the visualizations.
- Optional: input a name and location to save an image as the "outFilename" argument of the "heat_dendrogram" or "dendrogramDiff" functions.


Note: The length of the connections among nodes for dendrogram diagrams shown at the sides of the clustermaps, is given by:
- The distance measure d_ij between nodes, for the (Kt, CA, Antikt} algorithms 
- An integer starting at 0 and increasing by one each time a node is added, for the truth case. 

All the notebooks currently show visualizations of the same jet and the different options (truth, kt, CA, antikt)


## Jet visualizations comparison

![Fig.1](plots/1DTreeOnly/figTruth_jet10.png)

##### Fig. 1: ##### Fig. 1: Tree visualization of a sample jet generated with our model that represents a W boson jet. We show the values of **p** =(p_y,p_z) for each node and the scale Delta for the splitting of the inner nodes. The horizontal ordering of the leaves corresponds to the order in which the leaves are accessed when traversing the tree (and is not related to the particle momentum **p**).

Given a tree structure, this notebook gets heat clustermap visualizations of it. In particular, given a jet (the set of its constituents), visualization for the different trees obtained with each clustering algorithm (and the truth tree) are displayed. 

There are 2 functions for visualization:

a) function heat_dendrogram: Create  a heat dendrogram clustermap. If two jets are given as input:
- 1) If only one jet is given as input, both rows and columns are ordered according to that jet tree.
- 2) If truthJet and recluster_jet1, the visualization shows the heat data map of the truth jet, with columns ordered according to the truth jet and rows according to recluster_jet1.
- 3) If recluster_jet1 and recluster_jet2, the visualization shows the heat data map of recluster_jet2, with columns ordered according to the recluster_jet1 and rows according to recluster_jet2.

b) function dendrogramDiff:
Given two jet algorithms heat data matrices, we reorder the heat matrices according to the truth jet order (order in which leaves are accessed when traversing the truth tree) and take the difference.

Both functions have the same arguments:
- Truth jet dictionary
- recluster_jet1: reclustered jet 1
- recluster_jet2: reclustered jet 2
- full_path: Bool. If True, then use the total number of steps to connect a pair of leaves as the heat data. If False, then Given a pair of jet constituents {i,j} and the number of steps needed for each constituent to reach their closest common ancestor {Si,Sj}, the heat map scale represents the maximum number of steps, i.e. max{Si,Sj}.
- FigName: Dir and location to save a plot.


Usage:
- Input_dir: select the directory with the jets dictionaries.
- Input_jet: select the jet filename for the visualizations.
- Optional: input a name and location to save an image as the "outFilename" argument of the "heat_dendrogram" or "dendrogramDiff" functions.


Note: The length of the connections among nodes for dendrogram diagrams shown at the sides of the clustermaps, is given by:
- The distance measure d_ij between nodes, for the (Kt, CA, Antikt} algorithms 
- An integer starting at 0 and increasing by one each time a node is added, for the truth case. 


![Fig.1](plots/heatClustermap/figTruthTruth_singlepath_jet10.png)

##### Fig. 1: 


![Fig.2](plots/heatClustermap/figAntiktAntikt_singlepath_jet10.png)

##### Fig. 1:

![Fig.3](plots/heatClustermap/figDiffTruthKt_singlepath_jet10.png)

##### Fig. 1:




**Relevant Structure**:

<!--- [`data`](data/): Dir with the trees.-->
<!---->
<!---->
<!--<!---->-->
<!--<!--    -[`likelihood.py`](showerSim/likelihood.py): Calculate the log likelihood of a splitting node and of (a branch of) a tree. There are examples on how to run it at the end of the script.-->-->
<!---->
<!--- [`showerSim`](showerSim/): Dir with the simulation code.-->
<!---->
<!---[`exp2DShowerTree.py`](showerSim/exp2DShowerTree.py): Parton shower code to generate the trees. -->
<!---->
<!--- [`generate_jets`](scripts/generate_jets/):-->
<!---->
<!---[`generate_jets.py`](scripts/generate_jets/generate_jets.py): Calls and runs the parton shower code in [`showerSim`](showerSim/). The code could be run to get the augmented data as well.-->
<!---->
<!--<!---[`run2DShower.py`](showerSim/run2DShower.py): Run the parton shower code in [`showerSim`](showerSim/).-->-->
<!---->
<!--<!--- [`visualized-recursion_2D.ipynb`](visualized-recursion_2D.ipynb): Jet trees visualization.-->-->



##### **Running locally as a python package:**


<!--1. Clone the ToyJetsShower repository-->
<!--2. `cd ToyJetsShower`-->
<!--3. `make`-->



<pre>



</pre>

<img src="https://github.com/SebastianMacaluso/VisualizeBinaryTrees/blob/master/plots/IRIS-HEP.png" width="300" align="left"> <img src="https://github.com/SebastianMacaluso/VisualizeBinaryTrees/blob/master/plots/NYU.png" width="200" align="center"> <img src="https://github.com/SebastianMacaluso/VisualizeBinaryTrees/blob/master/plots/MSDSE.png" width="300" align="right">















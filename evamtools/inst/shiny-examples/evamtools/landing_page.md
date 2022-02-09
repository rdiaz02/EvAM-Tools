# What is Evam-tools?
***
```Evam-tools``` is a package that allows to infer cancer evolutionary pathways  from [_cross sectional data **(CSD)**_](#helpcsd) base on the output of [_cancer progression models  **(CPMs)**_](#cpms).

This web interface provides a user-friendly interactive version of the package. 
Here you can define complex scenarios with a few click and check the predictions of several models.

* In the ```Input``` tab you can define data sets to run ```Evam-tools```.
* In the ```Results``` tab you can see the output of [_CPMs_](#cpms). 

<center>
<img src="evamtools.png" width=850/>
</center>

```Evam-tools``` is is also available as an R package.
(See installation instructions in https://github.com/rdiaz02/EvAM-Tools#how-to-install-to-run-just-the-r-package)


# What is _cross-sectional data_?<a id="helpcsd"></a> 
***

In cross-sectional data, a single sample has been taken from each patient. That single sample represents the "observed genotype" of the tumor of that patient.
Genotype can refer to single point mutations, insertions, deletions or any other genetic modification.

In summary, _CSD_ is a binary matrix, filled either with _0_ if an event is not observed in a patient, or with _1_ if it is observed.

# How to use this web interface? <a id="input"></a>
***
The first step is to **define an scenario**.
You can do it by going to the ```User Input``` tab in the navigation bar at the top.

This web interface allows to define _CSD_ in three different ways:

* _By directly **defining genotypes frequencies**_ : you directly define what mutations are observed in how many patients.
* _By deriving genotype frequencies from a **direct acyclic graph (DAG)**_ : here you define dependency relationships between genes.
* _By deriving genotype frequencies from a **transition rate matrix**_ : the transition rate matrix reflects how genes affects each other, by making them more likely to mutate (positive theta) or less likely (negative theta).

Once you have created an scenario or selected one, you can hit the ```Run evamtools!``` button.
This will run several [_CPMs_](#cpms) and will display their [results](#helpresults). 

You can also increase or decrease the number of genes, or rename genes.

# How to build an scenario of cancer evolutions that makes some sense? A simple example <a id="examples"></a>

Image a simple scenarios where only study 3 genes (A, B ,and C). 

We sample 12 patients and we observe the following:
<center>

|            | A | B | C |
|:----------:|:-:|:-:|:-:|
| Patient  1 | 1 | 0 | 0 |
| Patient  2 | 1 | 0 | 0 |
| Patient  3 | 1 | 0 | 0 |
| Patient  4 | 1 | 0 | 0 |
| Patient  5 | 1 | 1 | 0 |
| Patient  6 | 1 | 1 | 0 |
| Patient  7 | 1 | 1 | 0 |
| Patient  8 | 1 | 1 | 0 |
| Patient  9 | 1 | 1 | 1 |
| Patient  10 | 1 | 1 | 1 |
| Patient  11 | 1 | 1 | 1 |
| Patient  12 | 1 | 1 | 1 |

</center>

Gene _A_ appers alone. Gene _B_ **always** appears when gene _A_ is mutated.
Hence, we can infeer thats the likelihood gene _B_ mutating is increased if gene _A_ is also mutated.
The same happens for gene _C_, it is only mutated if both genes _A_ and _B_ are also mutated.

In this example, we sample 16 patients and we observe the following:

<center>

|            | A | B | C |
|------------|---|---|---|
| Patient  1 | 1 | 0 | 0 |
| Patient  2 | 1 | 0 | 0 |
| Patient  3 | 1 | 0 | 0 |
| Patient  4 | 1 | 0 | 0 |
| Patient  5 | 0 | 1 | 0 |
| Patient  6 | 0 | 1 | 0 |
| Patient  7 | 0 | 1 | 0 |
| Patient  8 | 0 | 1 | 0 |
| Patient  9 | 1 | 1 | 0 |
| Patient  10 | 1 | 1 | 0 |
| Patient  11 | 1 | 1 | 0 |
| Patient  12 | 1 | 1 | 0 |
| Patient  13 | 1 | 1 | 1 |
| Patient  14 | 1 | 1 | 1 |
| Patient  15 | 1 | 1 | 1 |
| Patient  16 | 1 | 1 | 1 |

</center>

Genes _A_ and _B_ appears mutated either in combination or isolated. 
In principle, there is not a dependency relationship between those two

However, gene _C_ only appears mutated if _A_ and _B_ are also mutated.
We can infeer dependecy relationship:

>booth _A_ **AND** _B_ has to be mutated for mutation in gene _C_ to appear.

# How to interpret the ```Results```?<a id="helpresults"></a>
***

The results sections includes:

1. *Plotting the model*: here you can see either the **DAG** of each CPM showing the infered dependency relationships or the **transition rate matrix** in the case some [_CPMs_](#cpms).
1. *Plotting the sampling*: this plot is a bit more complex. It represents the flow that we have sampled using the output of its [_CPM_](#cpms). It highlights the most relevant genotypes and transitions between genotypes. It is created by making random samples using the parameters from the [_CPM_](#cpms) and counting the transitions observed and counting the transitions observed (defining edge width) and the genotype frequency (node size).  

<center>
<img src="transitions_dbxor.png" width=850>
</center>

3. *Tabular data*: represents the raw values computed from the model or extracted from the samples. This includes: 
  
  * *Transition rates*: the transition rate matrix of the continuous time Markov chain that models the transition from one genotype to another. This option is not available for OT, as OT does not return rates. <!-- or DBN. -->
  * *Transition probabilities*: conditional probability of transitions to a genotype (obtained using competing exponentials from the transition rate matrix; for OT this is actually a abuse of the untimed oncogenetic tree model, as explained in Diaz-Uriarte & Vasallo, 2019).
  * *Genotype frequencies*: absolute frequencies of each genotype as obtained by sampling from the given model. For all models except OT, obtained by simulating a sampling process from the transition rate matrix with observation time distributed as an exponential of rate 1. For OT, obtained from the code of Szabo & Boucher (package Oncotree) that gives the predicted probabilities of the genotypes according to the OT model; we then use multinomial sampling from the predicted probabilities. <!-- This option is not available for OT or DBN. -->

  * *Genotype transitions counts*: the number of times a transition from genotype A to genotype B has been observed when sampling. This option is not available for OT as this is undefined for OT. <!-- or DBN. -->

  * *Lambdas/probabilities*: parameters of each model. This option is not available for MHN, that returns hazards.
  * *Time-discretized transition matrix*:  the time-discretized version of the transition rate matrix. This option is not available for OT (as it requires rates).


# What is a cancer progression model (CPM)?<a id="cpms"></a>
***

Cancer progression models use cross-sectional data to infer probabilistic relationships between mutational events that lead to the disease. 

# What CPMs are included in ```Evam-tools```?<a id="cpms"></a>
***

*  **Oncogenetic Tress (OT):** this is the simplest graphical model. Restrictions are represented as a tree. Hence, a parent node can have many children, but children have a single parent.
*  **Conjuntive Bayesian Networks (CBN):** this model generalizes the tree-based restricion of OT to a direct acyclic graph (DAG). A DAG allows to include multiple parents. In CBN, when a node depends of many parent events, all of the them have to be present for the children to appear. In that sense, CBN models this relationships as conjuntive, in other words, it models the AND relationship.
<!-- *  **Disjuntive Bayesian Networks (DBN):** models multiple parent with the OR relationship. -->
<!-- *  **Monte Carlo CBN (MCCBN):** this is an implementation of CBN using Monte Carlo expectationi-maximization algorithm to work with a large number of mutations. -->
*  **Hidden Extended Suppes-Bayes Causal Networks (H-ESBCN):** is somewhat similar to CBN, but it includes automatic detection of logical formulas AND, OR and XOR. Note: H-ESBCN is used by its authors as part of Progression Models of Cancer Evolution (PMCE).
    
*  **Mutual Hazard networks (MHN):** in this model dpeendencies are not deterministic and events can make other events  more like or less likely (inhibiting influence). Hence, MHN includes multiple dependencies and is not limited to DAG schemes. The main parameters is a theta matrix of multiplicative hazards that represents how one event influences other events.

<center><h1>Evam-tools</h1></center>

&nbsp;
****
#### Table of contents <a id="evamtools"></a>
****
- [Evam-tools](#evamtools)
- [How to use this web interface?](#input)
<br>&nbsp;&nbsp;- [What is available under ```Results```?](#helpresults)
<br>&nbsp;&nbsp;- [Additional documentation](#additional_docs)
<br>&nbsp;&nbsp;- [Session timeouts, RAM and elapsed time execution limits](#limits)
<br>&nbsp;&nbsp;- [How long does it take to run](#timetorun)
- [What CPMs are included in ```Evam-tools```?](#cpms)
<br>&nbsp;&nbsp;- [Default options and default CPMs run](#cpms2run)
<br>&nbsp;&nbsp;- [References and related repositories](#refs)
- [Where is the code? Terms of use? Copyright](#code)
<br>&nbsp;&nbsp;- [Authors, contact and bug reports](#authors) 
- [Cookies](#cookies)


****
## Evam-tools <a id="evamtools"></a>
****

```Evam-tools``` is an R package and Shiny app that provides tools for evolutionary accumulation, or event accumulation, models. We use code from  "Cancer Progression Models" (CPM) but these are not limited to cancer (the key idea is that events are gained one by one, but not lost). ```Evam-tools``` is  also available as an R package (see https://github.com/rdiaz02/EvAM-Tools).


This web interface provides a user-friendly interactive version of the package. You can analyze your data, create cross-sectional data from scratch (by giving genotype frequencies), or generate data under different CPMs. You can compare results from different methods/models, as well as experiment and understand the consequences of changes in the input data on the returned inferences. You can also examine how a given model performs when data have been generated under other (or its own) model.

* In the ```User input``` tab (on top of the page) you can upload data or define cross-sectional data, or simulate cross-sectional data from models. These are then submitted to run.
* In the ```Results``` tab you can see the output.

<!-- You add/remove images by adding HTML code. The usual img block. But do not leave that commented code around or shiny might break. -->


&nbsp;
#### What is _cross-sectional_ data?<a id="helpcsd"></a> 

In cross-sectional data, a single sample is obtained from each subject or patient. That single sample represents the "observed genotype" of, for example, the tumor of that patient. Genotype can refer to single point mutations, insertions, deletions, or any other genetic modification. In this app, as is often done by CPM software, we store cross-sectional data in a matrix, where rows are patients or subjects, and columns are genes; the data is a 1 if the event was observed and 0 if it was not.

&nbsp;&nbsp;

****
## How to use this web interface? <a id="input"></a>
***


* Go first to the ```User input``` tab (on top of the page). Here you can:
    - Define your data by specifying the genotype composition or uploading a data set or a combination of both.
    - Generate data according to CPM models specified using DAGs (Directed Acyclic Graphs) and trees: CBN, OT, OncoBN, H-ESBCN: "DAG and rates/probs".
    - Generate data according to the MHN model: "MHN thetas".
    - When you generate data according to a model, you can specify the sample size, the amount of noise, if any, to add, and of course the parameters of the models.
    - You can also increase or decrease the number of genes, or rename genes.

* Change, if you want, options under "Advanced options and CPMs to use" (on the right of the screen).
* Click on "Run evamtools". 
* Results will be shown in the ```Results``` tab.


To make it easier to play with the tool, under "DAG and rates/probs" there are predefined DAGs that you can use to generate data; you can also modify the predefined DAGs before generating data.


&nbsp;&nbsp;

****
### What is available under ```Results```?<a id="helpresults"></a>
***

The results include:

* Plots:
    * On the top row, the DAGs and, for MHN, the log-&Theta; matrix.
        * The edges of the DAGs are annotated with the lambda (CBN, HESBCN),  weight (OT) or &theta; (OncoBN).
        * Remember: these are DAGs that have *genes (not genotypes) as nodes*. They represent the order restrictions of the events.
    * On the bottom row, the transition probabilities between genotypes or the transition rate matrix (for the models that return it); what is represented is chosen on the left, ```Customize the visualization```.
        * These plots have *genotypes (not genes) as nodes*.
        * (If you choose to *Sample for observed genotype transitions*, under  ```Advanced options and CPMs to use```, for models that return a transition rate matrix (CBN, H-ESBCN, MHN), you can also represent the observed genotype transitions.)
        * You can show, for the bottom plots, only some of the most relevant paths; again, modify options under ```Customize the visualization```.
        * The bottom plots might include genotypes never observed in the sample; these are shown in light green.
        * For easier visualization, in very busy plots, instead of the Genotypes you might want to show the last gene mutated.
    * You can represent only a subset of the fitted models (choose the CPMs to show). 
&nbsp;	
* Tabular output: a table where you can choose to display:
    * Transition probabilities: the probabilities of transition between genotypes. For OT and OncoBN see the additional documentation as these are not really available for untimed oncogenetic models.
    * Transition rates: for models that provide them (CBN, H-ESBCN, MHN) transition rates.
    * Predicted genotype relative frequencies: the predicted genotype frequencies from the fitted models.
    * Sampled genotype counts: the counts from obtaining a finite sample (of the size you chose) with the probabilities given by the predicted genotype frequencies. If you add noise, they include observational (e.g., genotyping) noise.
    * Observed genotype transitions (counts): if you choose to *Sample for observed genotype transitions* (under  ```Advanced options and CPMs to use```), for models that return a transition rate matrix (CBN, H-ESBCN, MHN), we obtain the observed sampled of genotypes by simulating sampling from the continuous time Markov chain; this provides also observed transition counts between genotypes.
&nbsp;		
* Original data: to help interpret the results, a histogram of the genotype counts is also provided.
&nbsp;
* You can also *download* the tabular results, fitted models, and the analyzed data.

&nbsp;&nbsp;
***
### Additional documentation<a id="additional_docs"></a>
***

Additional documentation is available from https://rdiaz02.github.io/EvAM-Tools/pdfs/Additional_doc_all.pdf.

(If you install the R package or the RStudio Docker image with the package, you also have access to the documentation of the package, which is included in this pdf).


&nbsp;&nbsp;



****
### Session timeouts, RAM and elapsed time execution limits<a id="limits"></a>
***

Inactive connections will timeout after 2 hours. The page will become gray, and if you refresh (e.g., F5 in most browsers) after this time, you will not get back your results, figures, etc, but start another session. 

Maximum RAM of any process is limited to 2 GB. Likewise, the analyses should be aborted after 1.5 hours of elapsed (not CPU ---we parallelize the runs) time. If you want to use the Shiny app without these limits, install a local copy. (To modify the time limit, change the value of variable EVAM_MAX_ELAPSED, in the definition of function "server", in file "server.R".  The RAM limit is imposed on the Docker containers we use; to remove it, run Docker without the memory limit.) Note: because of what we must do to enforce these limits, running over limits might not be signalled by an explicit error, but rather by a graying out or a complete refresh of the session.


&nbsp;&nbsp;

****
### How long does it take to run? <a id="timetorun"></a>
***

It depends on the number of features. For six features, and if you do not use H-ESBCN nor MC-CBN, it should take about 20 seconds.

&nbsp;&nbsp;


****
## What CPMs are included in ```Evam-tools```?<a id="cpms"></a>
***

*  **Oncogenetic Tress (OT):** Restrictions in the accumulation of mutations (or events) are represented as a tree. Hence, a parent node can have many children, but children have a single parent. OTs are untimed (edge weights represent conditional probabilities of observing a given mutation, when the sample is taken, given the parents are observed).

*  **Conjuntive Bayesian Networks (CBN):** This model generalizes the tree-based restriction of OT to a direct acyclic graph (DAG). A node can have multiple parents, and it denotes that all of the parents have to be present for the children to appear. Therefore, relationships are conjuntive (AND relationships between the parents). These are timed models, and the parameters of the models are rates given that all parents have been observed. We include both H-CBN as well as MC-CBN.

*  **Hidden Extended Suppes-Bayes Causal Networks (H-ESBCN):** Somewhat similar to CBN, but it includes automatic detection of logical formulas AND, OR and XOR. H-ESBCN is used by its authors as part of Progression Models of Cancer Evolution (PMCE). As for CBN, it returns rates.
  
*  **OncoBN**: Similar to OT, in the sense of being an untimed oncogenetic model, but allows both AND (the conjunctive or CBN model) and OR relationships (the disjunctive or DBN model).

  
*  **Mutual Hazard networks (MHN):** With MHN dependencies are not deterministic and events can make other events more like or less likely (inhibiting influence**. The fitted parameters are multiplicative hazards that represent how one event influences other events.

&nbsp;

****
### Default options and default CPMs run<a id="cpms2run"></a>
***

- In the Shiny app, by default we run CBN, OT, OncoBN, and MHN. If you want to run H-ESBCN or MC-CBN, or not run some of the above methods, (de)select them under ```Advanced options and CPMs to use```. (H-ESBCN or MC-CBN are not run by default, as they can take a long time).
- OncoBN can be run using a conjunctive or a disjunctive model. The default used in the Shiny app (and the ```evam``` function in the package) is the disjunctive model. You can use the conjunctive one by selecting it under ```Advanced options and CPMs to use```, in ```OncoBN options```, ```Model```.
- Most methods have other options that can be modified. Again, check ```Advanced options and CPMs to use```.





****
### References and related repositories<a id="refs"></a>
***

##### OT ####


- Desper, R., Jiang, F., Kallioniemi, O. P., Moch, H., Papadimitriou, C. H., &
  Sch\"affer, A A (1999). Inferring tree models for oncogenesis from comparative
  genome hybridization data. J Comput Biol, 6(1), 37–51.

- Szabo, A., & Boucher, K. M. (2008). Oncogenetic Trees. In W. Tan, & L. Hanin
  (Eds.), Handbook of Cancer Models with Applications (pp. 1–24). : World
  Scientific.

- Oncotree R package: https://CRAN.R-project.org/package=Oncotree

&nbsp;

##### H-CBN and MC-CBN ####

- Beerenwinkel, N., & Sullivant, S. (2009). Markov models for accumulating
  mutations. Biometrika, 96(3), 645.

- Gerstung, M., Baudis, M., Moch, H., & Beerenwinkel, N. (2009). Quantifying
  cancer progression with conjunctive Bayesian networks. Bioinformatics, 25(21),
  2809–2815. http://dx.doi.org/10.1093/bioinformatics/btp505


- Gerstung, M., Eriksson, N., Lin, J., Vogelstein, B., & Beerenwinkel, N. (2011). The Temporal Order of Genetic and Pathway Alterations in Tumorigenesis. PLoS ONE, 6(11), 27136. http://dx.doi.org/10.1371/journal.pone.0027136 

- Montazeri, H., Kuipers, J., Kouyos, R., B\"oni, J\"urg, Yerly, S., Klimkait, T., Aubert, V., … (2016). Large-scale inference of conjunctive Bayesian networks. Bioinformatics, 32(17), 727–735. http://dx.doi.org/10.1093/bioinformatics/btw459 
  
- GitHub repository for MC-CBN: https://github.com/cbg-ethz/MC-CBN

- Source code repository for h/ct-cbn:  https://bsse.ethz.ch/cbg/software/ct-cbn.html


&nbsp;

##### MHN ####


- Schill, R., Solbrig, S., Wettig, T., & Spang, R. (2020). Modelling cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1), 241–249. http://dx.doi.org/10.1093/bioinformatics/btz513 

- GitHub repository: https://github.com/RudiSchill/MHN  

&nbsp;

##### H-ESBCN (PMCE) ####


- Angaroni, F., Chen, K., Damiani, C., Caravagna, G., Graudenzi, A., & Ramazzotti, D. (2021). PMCE: efficient inference of expressive models of cancer evolution with high prognostic power. Bioinformatics, 38(3), 754–762. http://dx.doi.org/10.1093/bioinformatics/btab717 


-  Repositories and terminology: we will often refer to HESBCN, as that is the program we use, as shown here: https://github.com/danro9685/HESBCN. H-ESBCN is part of the PMCE procedure: https://github.com/BIMIB-DISCo/PMCE.

&nbsp;

##### OncoBN (DBN) ####

- Nicol, P. B., Coombes, K. R., Deaver, C., Chkrebtii, O., Paul, S., Toland, A. E., & Asiaee, A. (2021). Oncogenetic network estimation with disjunctive Bayesian networks. Computational and Systems Oncology, 1(2), 1027. http://dx.doi.org/10.1002/cso2.1027 

- GitHub repository: https://github.com/phillipnicol/OncoBN

&nbsp;

##### Conditional prediction of genotypes and probabilities of paths from CPMs ####

- Hosseini, S., Diaz-Uriarte, R., Markowetz, F., & Beerenwinkel, N. (2019). Estimating the predictability of cancer evolution. Bioinformatics, 35(14), 389–397. http://dx.doi.org/10.1093/bioinformatics/btz332 

- Diaz-Uriarte, R., & Vasallo, C. (2019). Every which way? On predicting tumor evolution using cancer progression models. PLOS Computational Biology, 15(8), 1007246. http://dx.doi.org/10.1371/journal.pcbi.1007246 

- Diaz-Colunga, J., & Diaz-Uriarte, R. (2021). Conditional prediction of consecutive tumor evolution using cancer progression models: What genotype comes next? PLOS Computational Biology, 17(12), 1009055. http://dx.doi.org/10.1371/journal.pcbi.1009055 

&nbsp;


****
## Where is the code? Terms of use? Copyright<a id="code"></a>
***

The complete source code for the package and the shiny app, as well information about how to run the shiny app locally, is available from https://github.com/rdiaz02/EvAM-Tools.


This app is free to use (if you have confidential data, you might want not to upload it here, and instead install the package locally). 


Most of the files for this app (and the package) are copyright Ramon Diaz-Uriarte and Pablo Herrera Nieto (and released under the Affero GPL v3 license ---https://www.gnu.org/licenses/agpl-3.0.html) except some files for HESBCN, MHN, and the CBN code; see full details in https://github.com/rdiaz02/EvAM-Tools#copyright-and-origin-of-files.


&nbsp;
#### Authors, contact and bug reports<a id="authors"></a>
Please, use the repository (https://github.com/rdiaz02/EvAM-Tools) and submit bug reports there. 


&nbsp;
****
## Cookies<a id="cookies"></a>
***
We use cookies to keep "sticky sessions" to the pool of servers (load balanced using [HAproxy](https://www.haproxy.org/)). By using the app, you confirm you are OK with this.



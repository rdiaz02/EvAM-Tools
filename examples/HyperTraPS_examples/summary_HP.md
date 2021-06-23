
## Set up
```r 
pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("hypertraps-process.R")
setwd(pwd0)
rm(pwd0)
```
---
## AND
### WT &#8594; A &#8594; B &#8594; C 

```r
dB_AND <- matrix(
  c(
    rep(c(1, 0, 0, 0), 100)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 1, 1, 0), 100)
    , rep(c(1, 1, 1, 1), 100)
    , rep(c(0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_AND) <- LETTERS[1:4]

do_HyperTraPS(dB_AND, "HP_AND", runs = 500, bi=200)
```

<img src="./HP_AND/freqs.jpg" height=300>
<img src="./HP_AND/stats-200.png" height=300>
<img src="./HP_AND/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_AND/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_AND/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_AND/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## OR
### WT &#8594; A &#8594; B &#8594; D 
### WT &#8594; A &#8594; C &#8594; D 

```r
dB_OR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200)
    , rep(c(1, 0, 1, 0), 100)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 1, 1, 0), 50)
    , rep(c(1, 1, 0, 1), 50)
    , rep(c(1, 0, 1, 1), 50)
    , rep(c(1, 1, 1, 1), 10)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_OR) <- LETTERS[1:4]

do_HyperTraPS(dB_OR, "HP_OR", runs = 500, bi=200)
```
<img src="./HP_OR/freqs.jpg" height=300>
<img src="./HP_OR/stats-200.png" height=300>
<img src="./HP_OR/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_OR/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_OR/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_OR/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## XOR
### WT &#8594; A &#8594; B &#8594; D
### XOR  
### WT &#8594; A &#8594; C &#8594; D

TODO: increase number of abcd
```r
dB_XOR <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200)
    , rep(c(1, 1, 0, 0), 100)
    , rep(c(1, 0, 1, 0), 100)
    , rep(c(1, 1, 0, 1), 50)
    , rep(c(1, 0, 1, 1), 50)
    , rep(c(0, 0, 0, 0), 10)
  ), ncol = 4, byrow = TRUE
)
colnames(dB_XOR) <- LETTERS[1:4]

do_HyperTraPS(dB_XOR, "HP_XOR", runs = 500, bi=200)
```
<img src="./HP_XOR/freqs.jpg" height=300>
<img src="./HP_XOR/stats-200.png" height=300>
<img src="./HP_XOR/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_XOR/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_XOR/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_XOR/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## OR then XOR
### A &#8594; B ; C &#8594; D; B XOR D for E 
Here I forgot to include some genotypes, like AC, but the algorithm is able to recover it. Hoewever, it fails also includes other genotypes that are not possible, like CE or ACE.

Check convergence plots
```r

dB_c1 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 0, 1, 1, 0), 100) #ACD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(1, 1, 1, 0, 1), 100) #ABCE
    , rep(c(1, 0, 1, 1, 1), 100) #ACDE
    , rep(c(1, 1, 1, 1, 0), 50) # ABCD
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c1) <- LETTERS[1:5]

do_HyperTraPS(dB_c1, "HP_c1", runs = 500, bi=200)
```
<img src="./HP_c1/freqs.jpg" height=300>
<img src="./HP_c1/stats-200.png" height=300>
<img src="./HP_c1/forwards-feature_graph_match-data_conditional_circular_graph.png" height=300>
<img src="./HP_c1/forwards-feature_graph_match-data_conditional_circular_adjacency.png" height=300>
<img src="./HP_c1/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c1/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## AND + OR
### ((A AND B) or C) to reach D
Here it aslo recovers missing genotypes like ACD or BCD.
```r 
dB_c2 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 0, 1), 200) #ABD
    , rep(c(0, 0, 1, 1), 250) #CD
    , rep(c(1, 1, 1, 1), 200) # ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c2) <- LETTERS[1:4]

do_HyperTraPS(dB_c2, "HP_c2", runs = 500, bi=200)
```
<img src="./HP_c2/freqs.jpg" height=300>
<img src="./HP_c2/stats-200.png" height=300>
<img src="./HP_c2/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_c2/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_c2/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c2/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## Separted paths
## WT &#8594; A AND B &#8594; E
## XOR 
## WT &#8594; C AND D &#8594; E
## ((A AND B) XOR (C AND D) to reach E

```r 
dB_c3 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 0, 0, 1), 100) #ABE
    , rep(c(0, 0, 1, 1, 1), 100) #CDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c3) <- LETTERS[1:5]

do_HyperTraPS(dB_c3, "HP_c3", runs = 500, bi=200)
```
<img src="./HP_c3/freqs.jpg" height=300>
<img src="./HP_c3/stats-200.png" height=300>
<img src="./HP_c3/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_c3/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_c3/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c3/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## Independent Paths
## WT &#8594; A &#8594; B
## WT &#8594; C &#8594; D
I forgot to include to include some genotypes,
but they are marked as probably missing ones.
```r 
dB_c4 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(0, 0, 1, 1), 200) #CD
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 0, 1, 1), 100) #ACD
    , rep(c(1, 1, 1, 1), 50) #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c4) <- LETTERS[1:4]

do_HyperTraPS(dB_c4, "HP_c4", runs = 500, bi=200)
```

<img src="./HP_c4/freqs.jpg" height=300>
<img src="./HP_c4/stats-200.png" height=300>
<img src="./HP_c4/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_c4/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_c4/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c4/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## Triple AND
## WT &#8594; (A AND B AND C) &#8594; D

```r 

dB_c5 <- matrix(
  c(
      rep(c(1, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0), 300) #C
    , rep(c(1, 1, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0), 200) #AC
    , rep(c(0, 1, 1, 0), 200) #BC
    , rep(c(1, 1, 1, 0), 100) #ABC
    , rep(c(1, 1, 1, 1), 50)  #ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c5) <- LETTERS[1:4]
do_HyperTraPS(dB_c5, "HP_c5", runs = 500, bi=200)

```
<img src="./HP_c5/freqs.jpg" height=300>
<img src="./HP_c5/stats-200.png" height=300>
<img src="./HP_c5/forwards-feature_graph_match-data_joint_circular_graph.png" height=300>
<img src="./HP_c5/forwards-feature_graph_match-data_joint_circular_adjacency.png" height=300>
<img src="./HP_c5/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c5/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

---
## Sign epistasis

Make another without A
```r 

## Sign epistasis
dB_se <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 300) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_se) <- LETTERS[1:2]
do_HyperTraPS(dB_se, "HyperTraPS_examples/HP_se", runs = 500, bi=200)
```

<img src="./HP_se/freqs.jpg" height=300>
<img src="./HP_se/stats-200.png" height=300>
<img src="./HP_se/forwards-feature_graph_match-data_conditional_circular_graph.png" height=300>
<img src="./HP_se/forwards-feature_graph_match-data_conditional_circular_adjacency.png" height=300>
<img src="./HP_se/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_se/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>
---

## Reciprocal Sign epistasis

```r 
dB_rse <- matrix(
  c(
      rep(c(1, 0), 50) #A
    , rep(c(0, 1), 50) #B
    , rep(c(1, 1), 400) #AB
    , rep(c(0, 0), 200) # WT
  ), ncol = 2, byrow = TRUE
)
colnames(dB_rse) <- LETTERS[1:2]
do_HyperTraPS(dB_rse, "HyperTraPS_examples/HP_rse", runs = 500, bi=200)
```

<img src="./HP_rse/freqs.jpg" height=300>
<img src="./HP_rse/stats-200.png" height=300>
<img src="./HP_rse/forwards-feature_graph_match-data_conditional_circular_graph.png" height=300>
<img src="./HP_rse/forwards-feature_graph_match-data_conditional_circular_adjacency.png" height=300>
<img src="./HP_rse/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_rse/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>
---

## ((A AND B) OR (C AND D) to reach E

```r 

dB_c7 <- matrix(
  c(
      rep(c(1, 0, 0, 0, 0), 300) #A
    , rep(c(0, 1, 0, 0, 0), 300) #B
    , rep(c(0, 0, 1, 0, 0), 300) #C
    , rep(c(0, 0, 0, 1, 0), 300) #D
    , rep(c(1, 1, 0, 0, 0), 200) #AB
    , rep(c(1, 0, 1, 0, 0), 200) #AC
    , rep(c(1, 0, 0, 1, 0), 200) #AD
    , rep(c(0, 1, 1, 0, 0), 200) #BC
    , rep(c(0, 1, 0, 1, 0), 200) #BD
    , rep(c(0, 0, 1, 1, 0), 200) #CD
    , rep(c(1, 1, 1, 0, 0), 100) #ABC
    , rep(c(1, 1, 0, 1, 0), 100) #ABD
    , rep(c(0, 1, 1, 1, 0), 100) #BCD
    , rep(c(1, 1, 0, 0, 1), 200) #ABE
    , rep(c(0, 0, 1, 1, 1), 200) #CDE
    , rep(c(1, 1, 1, 0, 1), 200) #ABCE
    , rep(c(1, 1, 0, 1, 1), 200) #ABDE
    , rep(c(1, 0, 1, 1, 1), 200) #ACDE
    , rep(c(0, 1, 1, 1, 1), 200) #BCDE
    , rep(c(1, 1, 1, 1, 1), 200) #ABCDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

do_HyperTraPS(dB_c7, "HP_c7", runs = 500, bi=200)
```

<img src="./HP_c7/freqs.jpg" height=300>
<img src="./HP_c7/stats-200.png" height=300>
<img src="./HP_c7/forwards-feature_graph_match-data_conditional_circular_graph.png" height=300>
<img src="./HP_c7/forwards-feature_graph_match-data_conditional_circular_adjacency.png" height=300>
<img src="./HP_c7/forwards-hypercube-graph-mach-data-g0.png" height=300>
<img src="./HP_c7/forwards-ws1-ws2.png" height=300>
<div style="page-break-after: always;"></div>

```r 
## Make a test with the real world data of the pmce paper

do_HyperTraPS("Bladder_Urothelial_Carcinoma.csv", "HyperTraPS_examples/HP_pmce", runs = 1000, bi=500)

WIP
```

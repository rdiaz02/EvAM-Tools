
En el plot de la derecha, la franja azul el walk simulator 1 (WS1): que busca ir desde 0 mutaciones a todo mutado. El naranja es el WS2, que busca solo busca dentro de los genotipos observados. PREGUNTA: Esto último es lo más útil para nuestros datos, no? 

```r 
pwd0 <- getwd()
setwd("../code_from_what_genotype_next/")
source("hypertraps-process.R")
setwd(pwd0)
rm(pwd0)

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
<img src="./HP_AND/stats-200.png" width=300>
<img src="./HP_AND/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_AND/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_AND/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_AND/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

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
<img src="./HP_OR/stats-200.png" width=300>
<img src="./HP_OR/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_OR/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_OR/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_OR/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

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
<img src="./HP_XOR/stats-200.png" width=300>
<img src="./HP_XOR/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_XOR/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_XOR/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_XOR/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>


```r
# A --> B ; C --> D; B XOR D for E 

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
<img src="./HP_c1/stats-200.png" width=300>
<img src="./HP_c1/forwards-feature_graph_match-data_conditional_circular_graph.png" width=300>
<img src="./HP_c1/forwards-feature_graph_match-data_conditional_circular_adjacency.png" width=300>
<img src="./HP_c1/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c1/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 

## ((A AND B) or C) to reach D

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
<img src="./HP_c2/stats-200.png" width=300>
<img src="./HP_c2/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_c2/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_c2/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c2/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
## WT --> A AND B --> E
## XOR 
## WT --> C AND D --> E
## ((A AND B) XOR (C AND D) to reach E
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
<img src="./HP_c3/stats-200.png" width=300>
<img src="./HP_c3/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_c3/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_c3/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c3/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
## WT --> A --> B
## WT --> C --> D
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

<img src="./HP_c4/stats-200.png" width=300>
<img src="./HP_c4/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_c4/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_c4/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c4/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
## WT --> (A AND B AND C) --> D

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

```
<img src="./HP_c5/stats-200.png" width=300>
<img src="./HP_c5/forwards-feature_graph_match-data_joint_circular_graph.png" width=300>
<img src="./HP_c5/forwards-feature_graph_match-data_joint_circular_adjacency.png" width=300>
<img src="./HP_c5/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c5/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
do_HyperTraPS(dB_c5, "HP_c5", runs = 500, bi=200)

## Sign epistasis
dB_c6 <- matrix(
  c(
      rep(c(1, 0, 0), 200) #A
    , rep(c(1, 1, 0), 20) #AB
    , rep(c(1, 1, 1), 400) #ABC
    , rep(c(0, 0, 0), 10) # WT
  ), ncol = 3, byrow = TRUE
)
colnames(dB_c6) <- LETTERS[1:3]
do_HyperTraPS(dB_c6, "HP_c6", runs = 500, bi=200)
```

<img src="./HP_c6/stats-200.png" width=300>
<img src="./HP_c6/forwards-feature_graph_match-data_conditional_circular_graph.png" width=300>
<img src="./HP_c6/forwards-feature_graph_match-data_conditional_circular_adjacency.png" width=300>
<img src="./HP_c6/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c6/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
## ((A AND B) OR (C AND D) to reach E

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
    , rep(c(1, 1, 0, 1, 1), 200) # ABDE
    , rep(c(1, 0, 1, 1, 1), 200) # ACDE
    , rep(c(0, 1, 1, 1, 1), 200) # BCDE
    , rep(c(1, 1, 1, 1, 1), 200) # ABCDE
    , rep(c(0, 0, 0, 0, 0), 10) # WT
  ), ncol = 5, byrow = TRUE
)
colnames(dB_c7) <- LETTERS[1:5]

do_HyperTraPS(dB_c7, "HP_c7", runs = 500, bi=200)
```

<img src="./HP_c7/stats-200.png" width=300>
<img src="./HP_c7/forwards-feature_graph_match-data_conditional_circular_graph.png" width=300>
<img src="./HP_c7/forwards-feature_graph_match-data_conditional_circular_adjacency.png" width=300>
<img src="./HP_c7/forwards-hypercube-graph-mach-data-g0.png" width=300>
<img src="./HP_c7/forwards-ws1-ws2.png" width=300>
<div style="page-break-after: always;"></div>

```r 
## Make a test with the real world data of the pmce paper

do_HyperTraPS("Bladder_Urothelial_Carcinoma.csv", "HyperTraPS_examples/HP_pmce", runs = 1000, bi=500)

WIP
```

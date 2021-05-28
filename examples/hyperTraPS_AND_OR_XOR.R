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

##Considerer here another examples to run

## ((A AND B) or C) to reach D

dB_c1 <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200) #A
    , rep(c(1, 1, 0, 0), 200) #B
    , rep(c(1, 1, 0, 0), 200) #C
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 0, 1, 0), 100)#AC
    , rep(c(1, 1, 0, 1), 50) #BC
    , rep(c(1, 0, 1, 1), 50) #ABC
    , rep(c(0, 0, 0, 0), 10) #ABD
    , rep(c(0, 0, 0, 0), 10) #CD
    , rep(c(0, 0, 0, 0), 10) # ABCD
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c1 <- LETTERS[1:4]

do_HyperTraPS(dB_c1, "HP_XOR", runs = 500, bi=200)
 
 

##Considerer here another examples to run

## ((A AND B) OR (C AND D) to reach E

dB_c1 <- matrix(
  c(
    rep(c(1, 0, 0, 0), 200) #A
    , rep(c(1, 1, 0, 0), 200) #B
    , rep(c(1, 1, 0, 0), 200) #C
    , rep(c(1, 1, 0, 0), 200) #E
    , rep(c(1, 1, 0, 0), 100) #AB
    , rep(c(1, 0, 1, 0), 100) #AC
    , rep(c(1, 1, 0, 1), 50) #AD
    , rep(c(1, 1, 0, 0), 200) #BC
    , rep(c(1, 1, 0, 0), 200) #BD
    , rep(c(1, 1, 0, 0), 200) #CD
    , rep(c(1, 0, 1, 1), 50) #ABC
    , rep(c(0, 0, 0, 0), 10) #ABD
    , rep(c(0, 0, 0, 0), 10) #BCD
    , rep(c(0, 0, 0, 0), 10) #ABE
    , rep(c(0, 0, 0, 0), 10) #CDE
    , rep(c(0, 0, 0, 0), 10) #ABCE
    , rep(c(0, 0, 0, 0), 10) # ABDE
    , rep(c(0, 0, 0, 0), 10) # ACDE
    , rep(c(0, 0, 0, 0), 10) # BCDE
    , rep(c(0, 0, 0, 0), 10) # ABCDE
    , rep(c(0, 0, 0, 0), 10) # WT
  ), ncol = 4, byrow = TRUE
)
colnames(dB_c1 <- LETTERS[1:4]

do_HyperTraPS(dB_c1, "HP_XOR", runs = 500, bi=200)

# A --> B ; C --> D

## ((A AND B) XOR (C AND D) to reach E
## Make a test with the real world data of the pmce paper

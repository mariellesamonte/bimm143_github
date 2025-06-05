# Class 6: R functions
Marielle Samonte (PID: A16861951)

- [1. Function basics](#1-function-basics)
- [2. Generate DNA sequence](#2-generate-dna-sequence)
- [3. Generate Protein Function](#3-generate-protein-function)

## 1. Function basics

Let’s start writing our first silly function to add some numbers:

Every R function has 3 things:

- name (we get to pick this)
- input arguments (there can be loads of these separated by a comma)
- the body (the R code that does the work)

``` r
add <- function(x, y=10, z=0){
  x + y + z
}
```

I can just use this function like any other function as long as R knows
about it (i.e. run the code chunk)

``` r
add(1, 100)
```

    [1] 101

``` r
add( x=c(1,2,3,4), y =100)
```

    [1] 101 102 103 104

``` r
add(1)
```

    [1] 11

Functions can have “required” input arguments and “optional” input
arguments. The optional arguments are defined with an equals default
value (`y=10`) in the function defination.

``` r
add(x=1,y=100,z=10)
```

    [1] 111

> Q. Write a function to return a DNA sequence of a user specified
> length. Call it `generate_dna()`

The `sample()` function can help here

``` r
#generate_dna <- function(size=5){ }

students <- c("jeff", "jeremy", "peter")

sample(students, size = 5, replace = TRUE)
```

    [1] "jeremy" "jeff"   "peter"  "jeff"   "peter" 

## 2. Generate DNA sequence

Now work with `bases` rather than `students`

``` r
bases <- c("A","C", "G", "T")
sample(bases, size = 10, replace = TRUE)
```

     [1] "C" "A" "A" "T" "C" "T" "A" "A" "T" "T"

Now I have a working ‘snippet’ of code T can use this as the body of my
first function version here:

``` r
generate_dna <- function(size=5) { 
  bases <- c("A","C", "G", "T")
  sample(bases, size = size, replace = TRUE)
  }
```

``` r
generate_dna()
```

    [1] "C" "C" "T" "T" "C"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(size=5, together=TRUE) { 
  bases <- c("A","C", "G", "T")
  sequence <- sample(bases, size = size, replace = TRUE)
  
  if(together) { 
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

``` r
generate_dna()
```

    [1] "GTGAT"

``` r
generate_dna(together = F)
```

    [1] "C" "C" "C" "C" "A"

## 3. Generate Protein Function

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length?

We can get the set of 20 natural amino-acids from the **bio3d** package.

``` r
aa <- bio3d::aa.table$aa1[1:20]
```

``` r
generate_protein <- function(size=6, together=TRUE) {
  
  ##Get the 20 amino acids as a vector
  aa <- bio3d::aa.table$aa1[1:20]
  sequence <- sample(aa, size = size, replace = TRUE)
  
  ## Optionally return a single element string
  if(together) {
    sequence <- paste(sequence, collapse = "")
  }
  return(sequence)
}
```

``` r
generate_protein()
```

    [1] "IQGQWK"

> Q. Generate random protein sequences of length 6 to 12 amino acids.

``` r
generate_protein(7)
```

    [1] "QPMAKVH"

``` r
generate_protein(8)
```

    [1] "WYHTAFIV"

``` r
generate_protein(9)
```

    [1] "KELIHMWYW"

``` r
# generate_protein(size=6:12)
```

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family of utility functions.

``` r
sapply(6:12, generate_protein)
```

    [1] "YKIKNI"       "VQPSEHN"      "QWFNRVST"     "ATRVVPNSG"    "CLREMPFDLI"  
    [6] "QWYNLRKMKNE"  "TMWGAKMRGAIL"

It would be cool and useful if I could get FASTA format output

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "LKFIFQ"       "SCKQKCR"      "QKQWFMAL"     "PVPDGAMKV"    "GCTVRDWVSL"  
    [6] "KNGFSMTVHRD"  "VWENIVNCNKPE"

``` r
cat(ans, sep="\n")
```

    LKFIFQ
    SCKQKCR
    QKQWFMAL
    PVPDGAMKV
    GCTVRDWVSL
    KNGFSMTVHRD
    VWENIVNCNKPE

I want this to look like FASTA format with an ID line. e.g.

    >ID.6
    YYKMTW
    >ID.7
    ERFAPFW
    >ID.8
    RYYCSTLL

The functions `paste()` and `cat()` can help us here…

``` r
cat(paste(">ID.", 6:12, "\n", ans, sep=""), sep="\n")
```

    >ID.6
    LKFIFQ
    >ID.7
    SCKQKCR
    >ID.8
    QKQWFMAL
    >ID.9
    PVPDGAMKV
    >ID.10
    GCTVRDWVSL
    >ID.11
    KNGFSMTVHRD
    >ID.12
    VWENIVNCNKPE

``` r
id.line <- paste(">ID.", 6:12, sep="")
id.line
```

    [1] ">ID.6"  ">ID.7"  ">ID.8"  ">ID.9"  ">ID.10" ">ID.11" ">ID.12"

``` r
id.line <- paste(">ID.", 6:12, sep="")
seq.line <- paste(id.line, ans, sep="\n")
cat(seq.line, sep="\n", file ="myseq.fa")
```

> Q. Determine if these sequences can be found in nature or are they
> unique

I BLASTp searched my FASTA format sequences against NR and found that
length 6, 7, 8 are not unique and can be found in the databases with
100% coverage and 100% identity.

Random sequences of length 9 and above are unique and can’t be found in
the databases.

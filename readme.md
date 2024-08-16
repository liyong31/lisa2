# Lisa2

Overview
=======

Lisa 2.0 is an LTLf to DFA conversion tool.
The previous version of Lisa features on a hybrid encoding for DFA construction, while this version only uses explicit state encoding. 

It is publicly available under the license GNU GPL v3.


Requirements
-----------------------------------

Lisa2 requires a C++14-compliant compiler.  G++ 5.x or later should work.

Third-party dependencies
-----------------------------------

* [Spot model checker version>=2.9](https://spot.lrde.epita.fr/)

* [MONA](https://github.com/liyong31/MONA)

Lisa2 is built on top of [Lisa](https://github.com/vardigroup/lisa) and relies on Spot and MONA to construct a DFA from a small LTLf formula.
When constructing a DFA from an LTLf formula with MONA, Lisa2 translates an LTLf formula to a formula in first order logic, which is then fed into MONA.

Complilation steps
=======

In the following we assume that we will compile Lisa on a Ubuntu system.

1. Install Spot

    Lisa needs Spot to convert an LTLf to a DFA and to perform intersection of DFAs with explicit state representation.

    * Get the latest version of Spot from https://spot.lrde.epita.fr/install.html.

    * Uncompress Spot and follow install intructions in README to install Spot. Note that the compilation of Spot may take a while.
    
            ./configure && make && sudo make install

    * Type ltl2tgba -f "F a" in command line, you expect to see some output starting with "HOA: v1".


2. Install MONA

    Lisa needs MONA to convert a formula in first order logic to a DFA.

    * Go to mona-1.4-17 directory and follow the install instructions in INSTALL.
    
            ./configure && make && sudo make install-strip

    **NOTE** 
    In BDD/bdd.h, the original table size is defined as #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000 (=2^24), which is too small for large DFA generation.
    We modify it to #define BDD_MAX_TOTAL_TABLE_SIZE 0x1000000000 (=2^36), so to allow MONA have larger table size during DFA construction.
    Note that MONA has explicit state representation but encodes the labels on transition symbolically.
    For more details on the representation of DFA in MONA, we refer to https://www.brics.dk/mona/mona14.pdf.
    
6. Compile Lisa2

    * Compile Lisa2 with Make:
    
            make T


Input format
=======

Lisa2 accepts LTLf formulas given as a .ltlf file written in SPOT format. 

For synthesis, it also requires a .part file. The .part file indicates the input and output propostitions for the synthesis task. 

Example .ltltf file

```
((COUNTER0 <-> INITCOUNTER0)) && (G (CARRY0 <-> INC)) && (G ((X COUNTER0 -> !(COUNTER0 <-> CARRY0)) && (X !COUNTER0 -> (COUNTER0 <-> CARRY0)))) && ((G ((!INC -> X INC)) -> F (!COUNTER0)))
```

Example .part file

```
.inputs INITCOUNTER0 INC
.outputs COUNTER0 CARRY0
```



For LTLf to DFA construction
==
To use the default setting to construct a DFA from an LTLf formula, type
    
        ./lisa2 -ltlf ./examples/ltlf3377.ltlf
    
    You are expected to see the output ending with "3377" for the number of states.



Syntax
==

The Linear Temporal Logic over finite traces (LTLf) has the same syntax as LTL.
Given a set P of propositions, the syntax of LTLf formulas supported by Spot is as follows:
```
φ ::= 1 | 0 | p | !φ | φ1 && φ2 | φ1 || φ2 | φ1 -> φ2 
      | φ1 <-> φ2 | X φ | X[!] φ | F φ | G φ | φ1 U φ2 | φ1 R φ2 | φ1 W φ2

```
where p ∈ P. Here 1 and 0 represent *true* and *false* respectively.
X (weak Next), X[!] (strong Next), F (Finally), G (Globally), U (Until), R (Release) and W (weak Until) are temporal operators.
We have that X[!] φ ≡ ! (X !φ), F φ = !(G !φ), φ1 U φ2 ≡ !(!φ1 R !φ2) and φ1 W φ2 ≡ G φ1 || (φ1 U φ2).
As usuall, we also have that F φ ≡ 1 U φ and G φ ≡ 0 R φ.

For the semantics of LTLf formula, we refer to [IJCAI13 paper](https://www.cs.rice.edu/~vardi/papers/ijcai13.pdf).
Specially, Spot supports a weak next and a strong next.
    
Weak next: *X a* is true if *a* holds at next step or if there is no next step.
In particular, *X(0)* is true iff there is no successor.
    
Strong next: *X[!] a* is true if *a* holds at next step and there must be a next step.
In particular *X[!]1* is true iff there is a successor.

**Note**
Other tools like [Syft](https://github.com/saffiepig/Syft) may interpret *X* as a strong next operator and use *N* to denote weak next operator.
Moreover, the minimal DFAs constructed by lisa may have one less state than those by MONA due to the fact that lisa removes nonaccepting sink state while MONA keeps it.

Please be noted the differences above.




Check for updates

## Multi-Grid Reaction-Diffusion Master Equation: Applications to Morphogen Gradient Modelling

Radek Erban1@ - Stefanie Winkelmann2

Received: 3 May 2024 / Accepted: 21 October 2024 / Published online: 27 November 2024 © The Author(s) 2024

## Abstract

The multi-grid reaction-diffusion master equation (mgRDME) provides a generalization of stochastic compartment-based reaction-diffusion modelling described by the standard reaction-diffusion master equation (RDME). By enabling different resolutions on lattices for biochemical species with different diffusion constants, the mgRDME approach improves both accuracy and efficiency of compartment-based reaction-diffusion simulations. The mgRDME framework is examined through its application to morphogen gradient formation in stochastic reaction-diffusion scenarios, using both an analytically tractable first-order reaction network and a model with a second-order reaction. The results obtained by the mgRDME modelling are compared with the standard RDME model and with the (more detailed) particle-based Brownian dynamics simulations. The dependence of error and numerical cost on the compartment sizes is defined and investigated through a multi-objective optimization problem.

Keywords Reaction-diffusion master equation · Multi-grid methods · Morphogen gradient formation · Stochastic simulation algorithms

## 1 Introduction

Compartment-based stochastic reaction-diffusion models have been used for modelling a range of biological processes, including Min-protein oscillations (Fange and Elf 2006) and ribosome biogenesis (Earnest et al. 2015) in E. coli, calcium signalling (Denizot et al. 2019; Dobramysl et al. 2016a), gene expression (Winkelmann and Schütte 2016), actin transport in filopodia (Zhuravlev et al. 2010; Erban et al.

| ख | Radek Erban                                                                                                     |
|---|-----------------------------------------------------------------------------------------------------------------|
|   | erban@maths.ox.ac.uk                                                                                            |
| ख | Stefanie Winkelmann                                                                                             |
|   | winkelmann@zib.de                                                                                               |
| 1 | Mathematical Institute, University of Oxford, Radcliffe Observatory Quarter, Woodstock Road, Oxford OX2 6GG, UK |
| 2 | Zuse Institute Berlin (ZIB), Takustrasse 7, 14195 Berlin, Germany                                               |

2014) and epidemic spreading (Winkelmann et al. 2021). They can be simulated using algorithms for continuous-time discrete-space Markov chains and are mathematically described using the reaction-diffusion master equation (RDME) for the probability mass function (Erban and Chapman 2020; Winkelmann and Schütte 2020). At the analytical level, the RDME builds a bridge between microscopic models for reaction-diffusion processes and macroscopic approximations through partial differential equations (Flegg et al. 2012; Montefusco et al. 2023; Flegg et al. 2015; Kang and Erban 2019).

To formulate a standard compartment-based model, the computational domain is discretized into compartments and diffusion is modelled as a jump process between the compartments (Winkelmann and Schütte 2020). Considering three-dimensional domains, it has been shown that the compartment size cannot be chosen arbitrarily small for systems containing second-order or higher-order reactions, i.e., the error of a compartment-based simulation increases as the compartment size approaches zero (Erban and Chapman 2020). There have been numerous studies (Isaacson 2009; Erban and Chapman 2009; Hellander et al. 2012; Kang et al. 2012; Isaacson 2013; Agbanusi and Isaacson 2014; Hellander and Petzold 2016; Isaacson and Zhang 2018) on the optimal choice of the compartment size and its influence on the approximation quality compared to microscopic modeling approaches given by the Doi model (Doi 1976a, b) or the Smoluchowski model (Andrews et al. 2010; Andrews 2012).

The optimal size of the compartment depends on the diffusion constant (Erban and Chapman 2009; Isaacson 2009; Hellander et al. 2012; Kang et al. 2012; Isaacson 2013; Agbanusi and Isaacson 2014; Hellander and Petzold 2016; Isaacson and Zhang 2018). In particular, if a biological system consists of molecular species with different diffusion constants, the compartment-based model can be naturally generalized to allow for different meshes (compartment sizes) for different chemical species (Cao and Erban 2014; Hellander and Hellander 2020; Li and Cao 2012). In this paper, we analyze such models using a generalization of the RDME which we will call the multi-grid reaction diffusion master equation (mgRDME). Numerical simulations of multi-grid reaction-diffusion models allow for high accuracy at reduced numerical cost compared to fully microscopically resolved systems (Cao and Erban 2014; Chatterjee et al. 2004; Dai et al. 2008). In some computational frameworks (Hellander and Hellander 2020), molecules can transfer from a fine-grained mesh to a coarse-grained mesh whenever appropriate. In the mgRDME that will be analyzed in this work, the grid size remains constant for each species but varies across species.

The mgRDME framework will be analyzed by applying it to biochemical systems of morphogen gradient formation (Saunders and Howard 2009; Kicheva et al. 2007; Bergmann et al. 2007; Gregor et al. 2007; Iyer et al. 2023). Morphogens are signaling molecules whose non-uniform distribution control pattern formation (or morphogenesis) during the development of multicellular organisms. Starting with the pioneering work of Turing (Turing 1952), it has been shown that pattern formation can naturally arise in reaction-diffusion systems through diffusion-driven instability, provided that the system contains at least two chemical species with different diffusion constants (Murray 2002). In particular, stochastic modelling of Turing patterns is a natural application area of the mgRDME framework (Cao and Erban 2014).

In this paper, we will focus on pattern formation (morphogen gradient formation) which results from pre-patterning, i.e., we assume that the studied domain has already been differentiated into two regions and the release of signalling molecules is localized in one of those regions (Rogers and Schier 2011; Shimmi and O’Connor 2003; Shimmi et al. 2005). The formation of morphogen gradients then leads to further pattern formation as the cells recognize and interpret the high or low morphogen concentration (Wolpert et al. 2015). Stochastic models helps us to understand the impact of noise (fluctuations) on the pattern (Teimouri and Kolomeisky 2018; Kolomeisky 2011; Erban and Chapman 2020).

In our model, we distinguish between signaling molecules A, which are produced locally in a specific spatial region, and morphogen molecules B, which are produced from A through a reaction that can occur anywhere in space. Both species move in space via diffusion. The locally restricted production of A (pre-patterning) leads to creation of a morphogen gradient which is then interpreted by cells to produce further patterning (Erban and Chapman 2020, Section 6.8). Since the signaling molecules  $A$  are assumed to diffuse at a significantly higher rate compared to  $B ,$  our reactiondiffusion systems are well suited to be analyzed via the mgRDME. Moreover, the primary focus lies on the spatial arrangement of the morphogen  $B ,$  represented by its gradient, while the specific spatial arrangement of signaling molecules A is of less importance. This prioritization is reflected in the use of distinct compartment sizes for the two species, with a finer resolution dedicated to the morphogen B.

We study two reaction networks. In Sect. 2, we consider a reaction-diffusion system with the two chemical species  $A$  (signal) and  $B$  (morphogen) being subject to the following three first-order chemical reactions

<span id="page-2-0"></span>

## ∅ −→ A (locally restricted) and A −→ B −→ ∅ , (1)

where the empty set ∅ is interpreted as sources and sinks of molecules (Erban and Chapman 2020), i.e., molecules of signal A are continuously produced and converted into molecules of morphogen  $B ,$  which are degraded over time. In Sect.3, we replace the first-order conversion reaction  $A \to B$  by the second-order  $d i m e r i z a t i o n$  reaction  $A + A \rightarrow B ,$  i.e., the chemical reactions (1) are changed to

<span id="page-2-1"></span>

## ∅ −→ A (locally restricted) and A + A −→ B −→ ∅, (2)

while the diffusion part of the reaction-diffusion model remains the same in both Sects. 2 and 3. Since morphogen gradient systems are effectively one-dimensional (with one ‘important’ direction), we study both systems in a one-dimensional domain. For reaction-diffusion systems with first-order chemical kinetics as in (1), the onedimensional results are directly applicable in higher dimensions. However, since our second system (2) includes the second-order reaction (dimerization), the diffusionlimited results will depend on the dimension of the physical space (Kang and Redner 1984; Montroll 1969; Ben-Avraham and Redner 1986).

For each of the two reaction networks (1) and (2), we formulate the standard RDME, the mgRDME, as well as a ‘ground truth’ model given by particle-based Brownian dynamics. The compartment-based modeling approaches are compared to

(a)

$$\operatorname { l e n g t h } L$$

$$( \mathrm { b } )$$

• A • B

0 0 0 6

tion

<span id="page-3-2"></span>Fig. 1 Standard compartment-based model. a Schematic of the computational domain [0, L] divided into  $K$  compartments of size  $h .$  In the first compartment on the left, highlighted in red, molecules of type  $A$  are produced at rate  $k _ { 1 } > 0 . \mathbf { b }$  The average number of molecules of  $A$  (red) and  $B$  (blue) at steady state obtained by solving Eq. (16) for parameters given in (18) and for  $L = 2$  and  $K = 4 0$  (i.e., we have  $h = L / K = 0 . 0 5 ,$   $\bar { D _ { A } } / h ^ { 2 } = 6 4$  and  $D _ { B } / h ^ { 2 } = 4 )$  (Color figure online)

the Brownian dynamics results based on the differences in their steady-state distributions. Specifically, we define an error function in terms of the distance between the distribution of morphogen molecules  $B$  at steady state, and a cost function as a measure of the numerical complexity when simulating the corresponding stochastic processes. Minimizing error and cost over possible choices of compartment sizes for  $A$  and  $B$  gives rise to a multi-objective optimization problem, which will be studied for both reaction networks. Our results are further summarized in the discussion Sect. 4. The time-dependent problems are discussed in Appendices 5.1 and 5.4.

<span id="page-3-0"></span>

## 2 First-Order Fast-Slow Morphogen Gradient Model

In this section, we consider a reaction-diffusion system with the first-order chemical reactions (1). We derive analytic expressions for the long-term spatial distribution of particles, examining the standard RDME, the mgRDME and the Brownian dynamics models in Sects.2.1, 2.2 and 2.3, respectively. Finally, we compare the models in Sect. 2.4 by means of a multi-objective optimization problem.

<span id="page-3-1"></span>

## 2.1 Standard Compartment-Based Model

We consider molecules of chemical species  $A$  and  $B$  that diffuse in the domain  $[ 0 , \, L \, \mathrm { l , }$  where  $L$  is the domain length. The diffusion coefficients are denoted  $D _ { A }$  and  $D _ { B } ,$ 

## Springer

respectively. They have physical units of [length]2/[time]. As for the standard RDME discretization, we divide the domain [0,  $L ]$  into  $K \ = \ L / h$  compartments of equal length  $h >$  0. This domain is schematically shown in Fig. 1a. Denoting the chemical species  $A$  (resp.  $B )$  in the  $i - t h$  compartment  $[ ( i \, \mathrm { ~ - ~ } 1 ) h , \, i \, h )$  by  $A _ { i }$  (resp.  $B _ { i } ) ,$  where  $i = 1 , \, 2 , \, \ldots , \, K ,$  then diffusion corresponds to two chains of ’chemicalreactions’’ (Cao and Erban 2014; Erban and Chapman 2020):

<span id="page-4-1"></span>

| −→←− dA \$A _ { 1 } \frac { d _ { 3 } } { d _ { A } } \, A _ { 2 } \, \frac { d _ { 3 } } { d _ { A } } \; . . . \frac { d _ { S } } { d _ { A } } \; A _ { K } , \qquad \qquad \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;\$ |  |  |     |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--|--|-----|
| −→ dB B1 −→←− d B B2 −→←− d B B3 ←− d ... ←− d BK ,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |  |  | (4) |

where dA = DA h2 and dB = DB h2 (5)

are the jump rates. Since  $h$  has the physical unit of [lengthl, the jump rates have units [time]−1. We assume that molecules of species A are produced in the first compartment on the left, highlighted in red in Fig. 1a. This corresponds to the reaction

<span id="page-4-0"></span>

## ∅ k1/h −→ A1 ,

(6)

where the rate constant  $k _ { 1 } \, > \, 0$  has units Itimel^{-1}. Scaling with  $h$  in reaction (6) is necessary to ensure that after multiplying by the compartment size  $h -$  as is usual for zero-order reactions (Erban and Chapman 2020) – we get a production rate independent of the compartment size. This gives the same production rate in the first compartment for any value of  $h ,$  and is consistent with our  $^ { \cdot }$  ground truth’ Brownian dynamics model corresponding to partial differential equations (PDEs) (27)–(28) in Sect. 2.3. Since (1) assumes that  $A$  is converted to  $B ,$  and  $B$  is degraded in the whole domain, we have the reactions  $A _ { i } \; \frac { k _ { 3 } } { \to } \; B _ { i } \; \; \; \; \; \mathrm { f o r } \; \; i = 1 , 2 , . . . . , K , \qquad ( \eta ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$ 

<span id="page-4-2"></span>with rate constants  $k _ { 2 }$  and  $k _ { 3 }$  again having physical units of  $[ \mathrm { t i m e l } ^ { - 1 } .$  The (random) number of particles in the  $i$  -th compartment at time  $t \geq 0$  is denoted by  $A _ { i } ( t )$  and  $B _ { i } ( t ) ,$  respectively. Let  $p ( \mathbf { n } , \mathbf { m } , t )$  be the joint probability that  $A _ { i } ( t ) = n _ { i } \, a n d \; B _ { i } ( t ) = m _ { i }$  for i = 1, 2,..., K, where we use the notation n =[n1, n2,..., nK ]∈ NK and  $\textbf { m } = [ m _ { 1 } , m _ { 2 } , \ldots , m _ { K } ] \in \mathbb { N } ^ { K }$  to define the system state. To formulate the RDME corresponding to the chemical reaction system (3)-(7), we define operators  ${ \cal O } _ { i } ^ { + } , { \cal O } _ { i } ^ { - }$  :  $\mathbb { N } ^ { K } \to \mathbb { N } ^ { K } \mathrm { ~ b y }$ 

$${ \cal O } _ { i } ^ { + } [ n _ { 1 } , \ldots , n _ { i - 1 } , n _ { i } , n _ { i + 1 } , \ldots , n _ { K } ] : = [ n _ { 1 } , \ldots , n _ { i - 1 } , n _ { i + 1 } , \ldots , n _ { K } ] , \quad ( 8 )
\mathcal { O } _ { i } ^ { - } [ n _ { 1 } , \ldots , n _ { i - 1 } , n _ { i } , n _ { i + 1 } , \ldots , n _ { K } ] : = [ n _ { 1 } , \ldots , n _ { i - 1 } , n _ { i + 1 } , \ldots , n _ { K } ] , \quad ( 9 )$$

for  $i \ = \ 1 , 2 , \ \ldots , \ K .$  By means of these operators, we define the diffusion operator  $\mathcal { D } : L ^ { 1 } ( \mathbb { N } ^ { K } \! \times \! \mathbb { N } ^ { K } ) \to L ^ { 1 } ( \mathbb { N } ^ { K } \! \times \! \mathbb { N } ^ { K } ) ,$  where

by

<span id="page-5-1"></span>
$$+ \frac { D _ { A } } { h ^ { 2 } } \sum _ { i = 1 } ^ { K } \left [ ( n _ { i } + 1 ) \, f ( { \cal O } _ { i } ^ { * } { \bf n } , { \bf m } ) - n _ { i } \, f ( { \bf n } , { \bf m } ) \right ]
+ \frac { D _ { B } } { h ^ { 2 } } \sum _ { i = 1 } ^ { \sum 1 } \left [ ( m _ { i } + 1 ) \, f ( \mathbf { n } , \mathcal { O } _ { i + 1 } ^ { + } \mathbf { m } ) - m _ { i } \, f ( \mathbf { n } , \mathbf { m } ) \right ]
\mathcal { R } f ( \mathbf { n } , \mathbf { m } ) : = k _ { 1 } \left [ f ( \mathcal { O } _ { 1 } ^ { - } \mathbf { n } , \mathbf { m } ) - f ( \mathbf { n } , \mathbf { m } ) \right ]$$

<span id="page-5-2"></span>and the reaction operator  $\mathcal { R } : L ^ { 1 } ( \mathbb { N } ^ { K } \times \mathbb { N } ^ { K } ) \to L ^ { 1 } ( \mathbb { N } ^ { K } \times \mathbb { N } ^ { K } )$  by

$$+ k _ { 3 } \sum \left \{ ( m _ { i } + 1 ) \, f ( { \bf n } , { \cal O } _ { i } ^ { + } { \bf m } ) - m _ { i } \, f ( { \bf n } , { \bf m } ) \right \} . \qquad \quad { ( 1 2 ) }$$

Then the reaction-diffusion master equation, which corresponds to the system of reactions (3-7), can be written as follows

{ . .1) = (0+) (s. m.). (13)

The stationary distribution is defined by

(n,m) ={imp(n,m, t)\$

and satisfies the stationary reaction-diffusion master equation

<span id="page-5-0"></span>0= (D + R) \phi(n, m). 14 (14)

Springer

## Solving (14), we obtain the product of Poisson distributions (Jahnke and Huisinga 2007)

<span id="page-6-3"></span>
$$\phi ( \mathbf { n } , \mathrm { m } ) = \exp \left [ - \sum _ { i = 1 } ^ { K } \left ( \bar { A } _ { i } + \overline { { B } _ { i } } \right ) \right ] _ { i = 1 } ^ { R } \frac { A _ { i } ^ { n } } { n _ { i } ! m _ { i } } , \qquad \qquad \qquad \qquad ( 1 5 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

where  $\overline { { A } } _ { i }$  and  $\overline { { B } } _ { i } ,$  for  $i = 1 , \, 2 , \ldots , \, K ,$  satisfy the steady-state equations

<span id="page-6-0"></span>

##  DA h2 S − k2 I  A =−k1e1,  DB h2 S − k3 I  B =−k2A. (16)

Here,  $I \in \mathbb { R } ^ { K \times K }$  is the identity matrix and vectors  $\overline { { \bf A } } \in \mathbb { R } ^ { K } , \overline { { \bf B } } \in \mathbb { R } ^ { K } , \textbf { e } _ { 1 } \in \mathbb { R } ^ { K }$  and matrix  $S \in \mathbb { R } ^ { K \times K }$  are given by

<span id="page-6-4"></span>

| A1                              | \$\overline { { B } } _ { 1 }\$    |       | −11 00 ... 00 ⎞ |                 |        |    |      |
|---------------------------------|------------------------------------|-------|-----------------|-----------------|--------|----|------|
| \$\overline { { A } } _ { 2 }\$ | \$\overline { { B } } _ { 2 }\$    | 0     | 1               | −210 ... 00     |        |    |      |
| A3 B                            | \$\overline { { B } } _ { 3 }\$ e1 | 0 S = | 0               | 1 −21 ... 00    |        |    |      |
| AK−1                            | BK−                                | 0     | 0               | 0               | 00 ... | −2 | 1    |
|                                 |                                    |       | 0               | 0 00 ... 1 −1 ⎠ |        |    |      |
|                                 |                                    |       |                 |                 |        |    | (17) |

The entries  $\overline { { A } } _ { i }$  and  $\overline { { B } } _ { i }$  are the average numbers of molecules of  $A$  and  $B ,$  respectively, in the i-th compartment at equilibrium. The solution (A,B) of Eq. (16) is plotted in Fig. 1b for the following choice of rate parameter values

<span id="page-6-1"></span> $k _ { 1 } = 1 0 0 , \quad k _ { 2 } = 2 , \qquad k _ { 3 } = 1 , \qquad D _ { A } = 0 . 1 6 , \qquad D _ { B } = 0 . 0 1 \, , \qquad ( 1 8 )$  showing the morphogen gradient. As the stationary distribution  $\phi$  is given by a product of Poisson distributions in (15), the values  $\overline { { A } } _ { i }$  and  $\overline { { B } } _ { i }$  not only determine the longterm averages but also the variances of the particle numbers at steady state. Since the values of  $\overline { { A } } _ { i }$  and  $\overline { { B } } _ { i }$  fully characterize the stationary distribution  $\phi ,$  we will compare the models in terms of  $\overline { { A } } _ { i }$  and  $\overline { { B } } _ { i } .$ 

<span id="page-6-2"></span>

## 2.2 Generalized Compartment-Based Model: mgRDME

Assuming that molecules of species A diffuse significantly faster than those of species  $B$  (as indicated by our parameter values in Eq. (18), where  $D _ { A } \gg \, D _ { B } ) ,$  we choose a larger compartment size for species A (Cao and Erban 2014), while keeping the compartments for species  $B$  of the same size as in the standard model in Sect.2.1. Let  $h _ { A } = L / K _ { A }$  for  $K _ { A } \in \mathbb { N }$  be the compartment length for species  $A ,$  and  $h _ { B } = L / K _ { B }$  for  $K _ { B } \in \mathbb { N } ,$  where  $K _ { B } \, > \, K _ { A } ,$  be the compartment length for species  $B .$  For simplicity, we will assume that the ratio  $\gamma \ : = \ K _ { B } / K _ { A } \ = \ h _ { A } / h _ { B } \ > \ 1$  is a natural number, i.e. γ ∈ N. However, a generalization is possible, see Remark 1 for an explanation. In the mgRDME formulation (Cao and Erban 2014), reaction (6) is replaced by

9 25 A.

(19)

which implies again that A-particles are produced at rate  $k _ { 1 }$  in the first compartment on the left (of the coarser A-discretization), independently of the grid size  $h _ { A } .$  The reaction  $A _ { i } \ \ \longrightarrow \ \ B _ { i }$  in (7) is generalized to

<span id="page-7-2"></span><span id="page-7-1"></span>A Bi, j = 1.2. Ka. for i () = (1, +1 ... j). (20) which means that after conversion from  $A$  to  $B$  the resulting  $B - m o l e c u l e i s p l a c e d$  uniformly in one of the y (smaller) compartments that overlap with the j-th (large) compartment of the reacting A-particle. The diffusion operator (11) is generalized to

while the reaction operator (12) now reads

<span id="page-7-0"></span>
$$\widetilde { \mathcal { R } } f ( \mathbf { n } , \mathbf { m } ) : = k _ { 1 } \left [ f ( \mathcal { O } _ { 1 } ^ { - } \mathbf { n } , \mathbf { m } ) - f ( \mathbf { n } , \mathbf { m } ) \right ]$$

$$+ \frac { k _ { 2 } } { \gamma } \sum _ { i = 1 } \sum _ { i = 1 } ^ { n } \left [ ( n _ { j } + 1 ) \, f ( { \cal O } _ { j } ^ { + } { \bf m } , { \bf m } ) - n _ { j } \, f ( { \bf n } , { \bf m } ) \right ]$$

Dividing by  $\gamma$  in the second line of (22) goes in line with summing over the set  $\mathcal { I } ( j )$  defined in (20) which contains  $| \mathcal { I } ( j ) | = \gamma$  elements; this can be interpreted as a sum of  $\gamma$  reactions  $A _ { j } \ \rightarrow \ B _ { i }$  that share the rate  $k _ { 2 } .$  In total, we obtain the  $m u l t i - g r i d$  reaction-diffusion master equation (mgRDME)

$$\ \frac { \partial p } { \partial t } ( \mathbf { n } , \mathbf { m } , t ) = ( \tilde { \mathcal { T } } + \tilde { \mathcal { R } } ) p ( \mathbf { n } , \mathbf { m } , t ) . \qquad \qquad \qquad ( 2 3 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

Solving the corresponding stationary mgRDME, given by  $0 = ( \bar { D } + \bar { \mathcal { R } } ) \, \phi ( \mathbf { n } , \mathbf { m } )$  in analogy to (14), we obtain a product of Poisson distributions similar to (15):

<span id="page-8-1"></span>where  $\overline { { A } } _ { j }$  and  $\overline { { B } } _ { i } ,$  for  $j \, = \, 1 , \, 2 , \, \ldots , \, K _ { A } , \,$  satisfy the generalized steady-state equations

 $\underline { { \underline { { H } } } } _ { A } \in \mathbb { R } ^ { K _ { A } \times K _ { A } }$  and  $I _ { B } \in \mathbb { R } ^ { K _ { B } \times K _ { B } }$  are identity matrices and vectors  $\overline { { \mathbf { A } } } \in \mathbb { R } ^ { K _ { A } } ,$  B {E}{b}, e1 RA and matrices SARAKA, Shad are defined in an analogous manner as before in (17). In addition, there isthe block matrix  $M \in \mathbb { R } ^ { K _ { B } \times K _ { A } }$  given by

 $M _ { 1 }$  100 0.. 0 1 0.. 0

with the non-zero entries (i.e., the ones) in the block  $M _ { j } ~ \in ~ \mathbb { R } ^ { \gamma \times K _ { A } }$  placed in the  $j$  -th column. As in the standard setting,  $\overline { { A } } _ { j }$  and  $\overline { { B } } _ { i }$  give both the long-term mean and variance of the population sizes in the respective boxes, and they fully characterize the stationary distribution  $\phi .$ 

Remark 1 The choice of an integer-valued ratio y simplifies notation in (20) and (22). However, a generalization to an arbitrary ratio  $K _ { B } / K _ { A } \ \not \in \ \mathbb { N }$  is straightforward: Let  $C _ { j } ^ { ( A ) } \, : = \left [ ( j - 1 ) h _ { A } , j h _ { A } )$  and  $C _ { i } ^ { ( B ) } \, : = \left [ ( i - 1 ) h _ { B } , i h _ { B } \right )$  denote the  $j - t h$  and  $i - t h$  compartment of the discretizations chosen for species  $A$  and species  $B ,$  respectively. Then, reaction  $A _ { j } \ \rightarrow \ B _ { i } ,$  meaning that a particle of species  $A$  located in  $C _ { j } ^ { ( A ) }$  is converted into  $B$  located in  $C _ { i } ^ { ( B ) } ,$  occurs at rate

(

<span id="page-8-0"></span>where |-| represents to the size of the compartment, specifically the length of the interval in the one-dimensional problems considered here. This means that the conversion rate  $k _ { 2 }$  is shared among the potential target compartments according to their overlap with the source compartment.

aled B

position  $\mathrm { p o s i t i o n }$  (a)  $\gamma = 2$ 

sition

(b\right)\gamma=4\$

<span id="page-9-0"></span>Fig. 2 Standard RDME versus mgRDME. Average number of molecules of A and B at equilibrium obtained by solving the standard stationary RDME (16) (red/blue dots) and by solving the stationary mgRDME (24) (grey bars) for rate constants given in (18) and for  $L = 2 , \, K _ { B } = 4 0$  and  $\mathbf { a } \, \gamma \, = \, 2 , \, \mathbf { b } \, \gamma = 4$  (Color figure online)

15 20

<span id="page-9-1"></span>Fig. 3 Difference between standard RDME and mgRDME. Difference  $||B1 – B3 ||1$  as a function of the ratio  $\, \gamma ~ = ~ K _ { B } / K _ { A } \,$  for rate constants given in (18) and for  $L ~ = ~ 2 , ~ K _ { B } ~ = ~ 6 0 ~ a n d ~ ~ K _ { A } ~ ~ \in ~$  [60, 30, 20, 15, 12, 10, 6, 4, 3], where  $B_{1$  the solution of the standard steady-state Eq. (16), and B\_{3, the solution of the generalized steady-state Eq. (24) (Color figure online)

## Model comparison: standard RDME versus mgRDME.

In the following, we denote the solution of the generalized steady-state Eq. (24) by  $\bar { \mathbf { A } } _ { \gamma } = ( \bar { A } _ { j } ^ { ( \gamma ) } ) _ { j = 1 , 2 , \ldots , K _ { A } } , \bar { \mathbf { B } } _ { \gamma } = ( \bar { B } _ { i } ^ { ( Y ) } ) _ { i = 1 , 2 , \ldots , K _ { B } }$  in order to emphasize its dependence on the ratio  $\gamma \, = \, K _ { B } / K _ { A }$  and to distinguish from the solution  $\bar { \bf A } = \bar { \bf A } _ { 1 } , \, \bar { \bf B } = \bar { \bf B } _ { 1 }$  of the standard compartment-based model given by Eq. (16). Figure 2 shows the solution  $\bar { \mathbf { A } } _ { \gamma } , \bar { \mathbf { B } } _ { \gamma }$  of equation (24) for the same parameter values as in Fig. 1b and for  $\gamma = 2$  and  $\gamma = 4 .$  Instead of  $\bar { \bf A } _ { \gamma }$  we plot the rescaled values  $\tilde { A } _ { i } ^ { ( \nu ) } : = \bar { A } _ { j } ^ { ( \nu ) } / \gamma$  for  $i \in \mathcal { I } ( j ) .$  This is for the purpose of comparability with the solution A of the standard steady-state Eq. (16) plotted in Fig. 1b.

We observe a qualitatively good agreement of the generalized solutions with the original one. To see how this agreement depends on the level of coarsening, we plot the difference  $\| \bar { \bf B } _ { 1 } - \bar { \bf B } _ { \gamma } \| _ { 1 }$  as a function of  $\gamma = K _ { B } / K _ { A }$  in Fig. 3. We choose  $\| \bar { \bf B } _ { 1 } - \bar { \bf B } _ { \gamma } \| _ { 1 }$  as the measure of error because the population of species  $B$  (morphogen) is the output of interest. Note that  $\bar { \bf B } _ { \gamma }$  is a vector of fixed size  $K _ { B }$  for any  $\gamma ,$  whereas the dimension of  $\bar { \bf A } _ { \gamma }$  is  $K _ { A } ,$  which depends on  $\gamma .$  Consequently, the distance between  $\bar { \bf A } _ { 1 }$  and  $\bar { \bf A } _ { \gamma }$  is not directly defined for  $\gamma \neq 1 .$  We see from Fig.3 that the difference  $\| \bar { \bf B } _ { 1 } - \bar { \bf B } _ { \gamma } \| _ { 1 }$  increases monotonically with γ , which is to be expected since larger values of γ represent a higher degree of coarse-graining. However, the increase of the error is nonlinear, and there is an area of values (γ - 5) where the error is sufficiently small to be deemed acceptable. Before we study the approximation quality in more detail, we define in the following section the ‘ground truth’ model given by particle-based Brownian dynamics.

<span id="page-10-0"></span>

## 2.3 ‘Ground Truth’ Particle-Based Model: Brownian Dynamics

As our ’ground truth’ model, we will consider a model with two species  $A$  and  $B$  diffusing and reacting in the unbounded domain R (Erban and Chapman 2020). Molecules of  $A$  are released at the origin,  $x \, = \, 0 , \,$  with rate  $2 k _ { 1 } ,$  which has units  $\operatorname { I s e c } ^ { - 1 }$  l. Due to symmetry, this is equivalent to studying the same model on the half-line [0,∞) with the rate constant  $k _ { 1 }$  and reflecting boundary conditions at  $x ~ = ~ 0$  (Erban and Chapman 2007). Molecules of  $A$  and  $B$  are diffusing with diffusion constants  $D _ { A }$  and DB, respectively. They are subject to the conversion and degradation reactions

| A k2 −→ B k3 −→ ∅. | (26) |
|--------------------|------|
|--------------------|------|

We compare the Brownian dynamics to the simulation on the finite interval [0, L], which is discretized into compartments in the RDME and mgRDME frameworks. The effect of the boundary at  $x = L$  is discussed in Sect.  $5 . 2$  of the Appendix.

Let  $a ( x , t )$  and  $b ( x , \, t )$  be the concentrations of molecules of  $A$  and  $B ,$  respectively, that is  $a ( x , \, t ) \, \mathrm { d } x$  denotes the average number of molecules in the interval  $[ x , \, x + \mathrm { d } x )$  at time  $t$  in the Brownian dynamics model. Then  $a ( x , \, t )$  and  $b ( x , \, t )$  satisfy the following PDEs (Kicheva et al. 2007, equation (1))

<span id="page-10-1"></span>
$$\frac { \partial a } { \partial t } = D _ { A } \frac { \partial ^ { 2 } a } { \partial x ^ { 2 } } - k _ { 2 } a + 2 k _ { 1 } \delta ( x ) , \qquad \qquad \qquad ( 2 7 ) \;$$

$$\frac { \partial b } { \partial t } = D _ { B } \frac { \partial ^ { 2 } b } { \partial x ^ { 2 } } + k _ { 2 } a - k _ { 3 } b , \qquad \qquad \qquad ( 2 8 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

where  $x \, \in \mathbb { R } ,$  and  $\delta ( x )$  is the Dirac delta function. We consider the boundary conditions

<span id="page-10-2"></span>

| lim a(x, t) = lim b(x, t) = 0. | (29) |
|--------------------------------|------|
| x→±∞ x→±∞                      |      |

Let  $\bar { a } ( x )$  and  $\bar { b } ( x )$  be the corresponding stationary distributions, i.e.

a¯(x) = lim t→∞ a(x, t), b¯(x) = lim t→∞ b(x, t).

They satisfy the steady-state equations

<span id="page-11-0"></span>0 = DA dx2 − k2 a¯ + 2k1 δ(x), (30) d2a¯ 0 = DB dx2 + k2 a¯ − k3 b¯ , (31) d2b¯

with boundary conditions

<span id="page-11-1"></span>
$$\operatorname* { l i m } _ { x \to \pm \infty } \bar { a } ( x ) = \operatorname* { l i m } _ { x \to \infty } \bar { b } ( x ) = 0 . \qquad \qquad \qquad ( 3 2 )$$

We notice that the PDE (27) (resp. the ODE (30)) does not depend on  $b ( x , t )$  (resp. b(x)). In what follows, we will therefore first solve the Eq. (27) for the chemical species A, and then we will substitute it into the second Eq. (28) for the chemical species  $B$  to get our  $^ { \cdot }$  ground truth’ solution. In the following, we will do this for the steady-state problem. The time-dependent problem is considered in the Appendix, see Sect. 5.1.

Steady-state solution

Since the second derivative of the absolute value function, |x|, is equal to 2δ(x), we can write the solution to Eq. (30) with boundary conditions (32) as

<span id="page-11-3"></span>
$$\tilde { a } ( x ) = \frac { k _ { 1 } } { \sqrt { D _ { A } \, k _ { 2 } } } \exp \left [ - \sqrt { \frac { k _ { 2 } } { D _ { A } } } \, | \, , \qquad \qquad \qquad \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

This result is also obtained when taking the limit  $t \ \rightarrow \ \infty$  of the time-dependent solution  $a ( x , t )$  as given in Eq. (49) of Appendix 5.1. Substituting into Eq. (31), we get

DB dx2 + k1 d2b¯ DA k2 exp − DA k2 |x| − k3 b¯ = 0.

Consequently, using boundary conditions (32), we obtain

<span id="page-11-2"></span>b¯(x) = k1 DA k3 − DB k2 √k2 DA exp − DA k2 |x| + (DB k2 − DA k3) k1k2 √DB √k3 exp − DB k3 |x| (34)

for the concentration  $\bar { b }$  of molecules of  $B \,$  at equilibrium. Our compartment-based models simulate the problem on the finite interval [0, L], with zero-flux boundary condition at  $x = L .$  To get the  $^ { \cdot }$  ground truth’ solution (34), we have used the infinite domain, R, with boundary conditions (32). In Sect. 5.2 we show that the effect of non-flux boundary conditions on the  $^ { \mathrm { ~ s ~ } }$  ground truth’ solution  $\bar { a }$  vanishes for  $L \rightarrow \infty .$  In particular, the concentration is negligible at  $x \, \geq \, L$  as long as  $L$  is large enough, and we can use (33) and (34) for a comparison with the compartment-based models.

<span id="page-12-0"></span>

## 2.4 Model Comparison

Since the reaction system under consideration is a first-order reaction network, the standard compartment-based model of Sect. 2.1 converges to the Brownian dynamics solution for  $h \ \rightarrow \ 0$  (Engblom et al. 2009; Erban and Chapman 2020). Our goal is now to compare Brownian dynamics with the compartment-based model for finite h, both for the standard RDME and for the mgRDME. In application, one is mainly interested in the output given by the product  $B ,$  which motivates to define the distance between the models by means of the average number of molecules of B.

## 2.4.1 Error Definition

Given the steady-state solution  $\bar { \bf B } = ( \bar { B } _ { 1 } , \, \bar { B _ { 2 } } , \ldots , \, \bar { B } _ { K _ { B } } )$  of the mgRDME model with  $K _ { B }$  boxes of size  $h = L / K _ { B }$  (which agrees with the standard RDME solution when choosing  $\gamma = 1 ) ,$  we define the piecewise constant function  $\bar { b } _ { h } :$  [0,  $L ] \to$  I0,  $\infty ) ,$ 

<span id="page-12-2"></span>
$$\bar { b } _ { h } ( x ) : = \bar { B } _ { i } / h \quad \mathrm { f o r } \, x \in [ ( i - 1 ) h , i h ) , \; \; i = 1 , 2 , \ldots , K _ { B } , \qquad ( 3 5 )$$

which approximates the steady-state concentration of molecules of  $B$  as it depends on the location  $x .$  The  $L ^ { 1 }$  -distance to the Brownian dynamics solution  $\bar { b }$  defines

<span id="page-12-1"></span>K err(KB,γ) := ||b¯ − b¯ h||L1 = 0 L |b¯(x) − b¯h(x)| dx =  i=1 B (i−1)h ih |b¯(x) − B¯i| dx , (36)

the spatial error in dependence on  $K _ { B }$  and on the ratio  $\gamma \equiv K _ { B } / K _ { A } .$  We note that this error would be positive even if  $\bar { b } _ { h }$  was a discretization of  $\bar { b }$  (and not the solution of the compartment-based model), simply because of the coarse-graining. We thus expect the error to naturally decrease with decreasing h.

Figure 4 shows this error in dependence on the grid size  $h ,$  where  $h _ { A } = h _ { B } = h$  for the standard RDME, while  $h _ { A } = h$  and fixed  $h _ { B } = h ^ { * } = 1 / 6 0$  for the generalized mgRDME. For both scenarios, there is a monotone increase of the error with h, which is to be expected in a first-order reaction system. The generalized model shows better results as compared to the standard one, which is due to the fixed small grid size  $h ^ { * }$  for  $B .$  However, the error  $\operatorname { e r r } ( K _ { B } , \, \gamma )$  seems to be linear in  $h$  for  $\gamma = 1$  (standard model, blue crosses), while it is nonlinear in  $h = h _ { A }$  for the generalized system and stays close to zero for small  $h = h _ { A } .$  We can say that for  $h = h _ { A } \leq 0 . 1 \ ( \gamma \ \leq 6 )$  the generalized

| rd        |     |     |
|-----------|-----|-----|
| eralized  |     |     |
| 05        |     |     |
| 0.2 \$h\$ | 0.3 | 0.4 |

<span id="page-13-0"></span>Fig. 4 Error in dependence on the compartment size. Error err(K g, y), defined in (36), between the steadystate solution  $\bar { b }$  of the  $^ { \cdot }$  ground truth’ model given by (34) and its approximation by the standard and the generalized compartment based models with solutions  $\bar { b } _ { h }$  defined in (35) in dependence on the compartment size h. We have h 4 = h g = h\$ for the standard model, while h 4 = h and h g = h^{k} = 1/60 for the multi-grid model. Rate constants given in (18) and  $L = 2$  (Color figure online)

mgRDME model gives a good approximation, while the standard compartment-based model already shows a clear error.

Despite the availability of an analytical solution for the steady-state distribution, there remains a fundamental interest in stochastic simulations of the dynamics, which allow for a more comprehensive exploration of the system’s behavior including complex interactions and rare events. We point out that, in Fig.4 and for a given  $h ,$  the numerical effort for simulating trajectories is larger for the multi-grid model compared to the standard model due to the retention of the small grid size  $h _ { B } = h ^ { * }$  for  $B .$  However, the overall goal of the generalized method is to reduce numerical cost while simultaneously minimizing errors. Therefore, determining an optimal value for y entails multi-objective optimization, a task we will address next.

## 2.4.2 Multi-Objective Optimization

In this section, we study the dependence of the spatial error err(Kp, y) defined in (36) and of the numerical effort on the grid sizes  $h _ { A }$  and  $h _ { B } .$  We calculate error and cost for different parameter values to show the advantage of the generalized method over the standard method. To define the cost function, we use propensities at steady state because these indicate the runtime of our numerical simulations.

Propensities at steady state

The total number of particles at steady state is given by,1

The overall propensity for a reaction to take place at steady state (no matter which reaction or where in space) is

<span id="page-13-2"></span>
$$\bar { r } _ { \mathrm { r e a c } } : = k _ { 1 } + k _ { 2 } \, \bar { A } _ { \mathrm { t o t a l } } = 3 \, k _ { 1 } \, . \qquad \qquad \qquad ( 3 8 ) \, \,$$

<span id="page-13-1"></span>

 $^ { 1 }$  The values resultfrom setting the related reactionrate equations to zero, i.e. from solving  $k _ { 1 } - k _ { 2 } \bar { A } _ { \mathrm { t o t a l } } = 0$  and  $k _ { 2 } \bar { A } _ { \mathrm { t o t a l } } - k _ { 3 } \bar { B } _ { \mathrm { t o t a l } } = 0 .$  We note that this approach is only applicable to first-order reaction systems.

Analogously, the overall propensity for a diffusive jump to take place at steady state is

r¯diff := h2 DA (2A¯total − A¯1 − A¯KA ) + h2 DB (2B¯total − B¯1 − B¯KB )

<span id="page-14-1"></span>- 2 h2 DA A A¯total + 2 h2 DB B B¯total = 2 h2 DA A k2 k1 + 2 h2 DB B k1 k3 (39)

with the negative terms in the first line resulting from the fact that diffusive jumps at the outer boundaries of the domain can only go in one direction. The latter expression is a good approximation as long as the population sizes in the outer boxes are not too large (as compared to the total population size), which is the case in the scenarios considered.

The number of iterations steps needed to create a trajectory of the dynamics (starting in steady state, using the Gillespie algorithm (Gillespie 1977, 2007)) scales with  $\bar { r } _ { \mathrm { r e a c } } + \bar { r } _ { \mathrm { d i f f } } .$  This motivates to define the cost function as a sum of (38) and (39), giving

<span id="page-14-2"></span>
$$c ( K _ { B } , \gamma ) : = 3 k _ { 1 } + 2 \frac { D _ { A } } { h _ { A } \, k _ { 2 } } \frac { k _ { 1 } } { k _ { 3 } } , \qquad \qquad \qquad ( 4 0 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

where  $h _ { B } ~ = ~ L / K _ { B } ~$  and  $h _ { A } ~ = ~ \gamma h _ { B } .$  Figure5 shows the values of error and cost depending on  $K _ { B }$  and  $\gamma$  (where  $\gamma \, = \, 1$  corresponds to the standard RDME model). We see that the generalized model clearly outperforms the standard model: The Pareto front exclusively comprises points derived from the generalized model. Consequently, the decrease in spatial resolution for species A introduced by the generalized model minimally impacts the error in the  $B - s o l u t i o n , r e l a t i v e t o t h e a l t e r a t i o n i n n u m e r i c a l$  cost resulting from this coarse-graining. We note the special role of the value  $\gamma = 4 =$   $\sqrt { D _ { A } / D _ { B } }$  which can be considered a good starting point for a practicalimplementation of the mgRDME as it balances the jump rates of the two species and enhances the model’s efficiency.

Remark 2 We note that the scenario would change if we were to define an error metric that incorporates the distance within the A-population. This is evident because ||A −  $\bar { \mathbf { A } } ( K _ { B } , \, \gamma ) | | _ { 1 }$  rises with  $\gamma \, ,$  irrespective of  $K _ { B } .$ 

In the examined first-order reaction system, only the production of A depends on location, whereas the other reactions occur independently of the molecules’ spatial position. A stronger effect of the grid size on the reaction dynamics is expected for a second-order reaction system, which is studied next.

<span id="page-14-0"></span>

## 3 Fast-Slow Morphogen Gradient with Dimerization

We consider the fast-slow morphogen gradient system of Sect. 2, replacing the conversion reaction  $A \rightarrow B$  by the second-order reaction of  $d i m e r i z a t i o n \ A + A \to B ,$  i.e., the chemical reactions (1) are changed to (2), while the diffusion part of the reactiondiffusion model remains consistent with that described in the previous section. Again,

+ standard

| +                            | standard |
|------------------------------|----------|
| generalized with K_{0}= 120  |          |
| O generalized with K_{0j}=60 |          |

|      | od 11   |     |       |   |      |       |    |      |     |
|------|---------|-----|-------|---|------|-------|----|------|-----|
|      | \$60 11 |     |       |   |      |       |    |      |     |
| 1202 |         |     |       |   |      |       |    |      |     |
|      |         | 0 1 |       |   |      |       |    |      |     |
|      |         |     |       |   |      |       |    |      |     |
| § 10 |         |     |       |   |      |       |    |      |     |
|      |         |     |       |   | 10 1 |       |    |      |     |
|      |         |     | 10 21 |   |      |       |    |      |     |
|      |         |     |       |   |      |       |    | 0 11 |     |
|      |         |     | F603) |   |      |       |    |      |     |
|      |         |     | 1     | 1 | 804  | 1     | on |      |     |
| 1    | 1.5     | 2   | 2.5   | 3 |      | 3.5   | 4  | 4.5  | 5.5 |
|      |         |     |       |   |      | error |    |      |     |

<span id="page-15-0"></span>

Fig. 5 Error and cost depending on the method. Error err(Kg, y) as defined in (36) and numerical cost  $c ( K _ { B } , \gamma )$  as defined in (40) for different values of  $K _ { B }$  and  $\gamma \, = \, K _ { B } / K _ { A }$  (given by labels of the form [K\_}\$ y]). The colored symbols are for orientation: blue pluses for the standard RDME, red stars for the generalized mgRDME with  $K _ { B } = 1 2 0$  and green circles for the generalized mgRDME with  $K _ { B } = 6 0 .$  Parameters are given in (18) and  $L = 2$  (Color figure online)

we choose the one-dimensional domain [0, L]. The associated modeling approaches are introduced in Sect. 3.1 and compared in Sect. 3.2.

<span id="page-15-1"></span>

## 3.1 Modeling Approaches for the Dimerization System

Considering the compartment-based models, the diffusion operators are the same as in Sect.2, i.e., for the standard RDME model it is giyen by  $\mathcal { D }$  as defined in (11), and for the generalized mgRDME model it is given by  $\bar { \mathcal D }$  as defined in (21). In contrast, the reaction operators now explicitly depend on the grid size  $h$  via the second-order reaction  $A + A \to B .$  For the standard compartment-based model, given by the RDME  $\frac { \partial p } { \partial t } ( \mathbf { n } , \mathbf { m } , t ) = ( \mathcal { D } + \mathcal { R } ) \ p ( \mathbf { n } , \ \mathbf { m } , \ t ) ,$  it reads

$$\mathcal { R } f ( \mathbf { n } , \mathbf { m } ) : = \, k _ { 1 } \left [ f ( \mathcal { O } _ { 1 } ^ { - } \mathbf { n } , \mathbf { m } ) - f ( \mathbf { n } , \mathbf { m } ) \right ]$$

$$+ \frac { k _ { 2 } } { h } \sum \left [ ( n _ { i } + 2 ) ( n _ { i } + 1 ) \, f ( { \cal O } _ { i } ^ { + } { \bf \rho } _ { i } ^ { - } { \bf m } ) - \operatorname* { m a x } \{ n _ { i } ( n _ { i } ( n _ { i } - 1 ) , 0 \} \, f ( { \bf n } , { \bf m } ) \right ]$$

$$+ k _ { 3 } \sum _ { i = 1 } ^ { K } \left [ ( m _ { i } + 1 ) \, f ( \mathbf { n } , \mathcal { O } _ { i } ^ { + } \mathbf { m } ) - m _ { i } \, f ( \mathbf { n } , \mathbf { m } ) \right ]$$

for  $h \, = \, h _ { A } \, = \, h _ { B } .$  The rate constant  $k _ { 2 }$  now has physical units Ilengthl/timel, i.e., units of the  $^ { \cdot }$  one-dimensional volume’ divided by time, because  $k _ { 2 }$  is the rate constant of a bimolecular reaction, in contrast to  $k _ { 2 }$  in (12), which has physical units Itimel^{-1}. Dividing by  $h$  in the second line stems from the standard scaling of a second-order reaction rate by the volume of the domain in which it takes place, which is given by  $h$  for the one-dimensional case under consideration. Similarly, for the generalized mgRDME model, we have

$$\tilde { \mathcal { R } } f ( \mathbf { n } , \mathbf { m } ) \ : = \ k _ { 1 } \left [ f ( \mathcal { O } _ { 1 } ^ { - } \mathbf { n } , \mathbf { m } ) - f ( \mathbf { n } , \mathbf { m } ) \right ]$$

$$+ k _ { 3 } \sum _ { i = 1 } ^ { K _ { B } } \left [ ( m _ { i } + 1 ) \, f ( { \bf n } , { \cal O } _ { i } ^ { + } { \bf m } ) - m _ { i } \, f ( { \bf n } , { \bf m } ) \right ] .$$

The factor  $1 / \gamma$  in the second summand results from splitting the reaction  $A + A \to$   $B$  occurring in the  $j - t h$  compartment into  $| \mathcal { I } ( j ) | = \gamma$  possible ones, depending on the placement of the product  $B \colon$ 

## Aj + y 2 Bi, j = 1.2. Ka. i () = (1, +1, ... jy), (41)

## which conforms with the first-order case, see (20) and (22).

## ’Ground truth’ particle-based model: Brownian dynamics

As the ’ground truth’we choose particle-based dynamics given by the Doi or \lambda-0 model (Doi 1976a, b; Erban and Chapman 2009; Lipkova et al. 2011). Particles move in space by Brownian motion, and two particles of species A undergo dimerization  $A + A \rightarrow B$  at rate  $\lambda \, > \, 0$  whenever being within a separation  $\varrho \, > \, 0 ,$  called the reaction radius.

## Relation between rate constants for bimolecular reactions

Let  $h$  denote the grid size of the compartment-based model (given by the standard RDME), and let  $k _ { 2 } / h ^ { d }$  (where  $d$  is the dimension of physical space) denote the rate for two molecules of A to react when located in the same compartment. The (standard) RDME may be seen as a formal approximation of the Doi model (Isaacson 2008), but it loses second-order reactions as \$h \ \rightarrow \ 0 \ \$ (Isaacson and Peskin 2006; Isaacson 2009; Erban and Chapman 2009). In Isaacson (2013) a convergent RDME has been developed, which allows second-order reactions of molecules located in different compartments and thereby does not lose bimolecular reactions as  $h \rightarrow$  0. The standard RDME may be interpreted as an asymptotic approximation of the convergent RDME for  $\varrho / h \ll 1$  (Isaacson 2013), where  $\varrho$  is the reaction radius of the Doi model. In this case of a comparatively large compartment size, particles can be expected to only react when being located in the same box, and an adequate choice of the binding rate constant  $k _ { 2 }$  is given by  $k _ { 2 } \, = \, \lambda | B _ { \circ } | ,$  where  $| B _ { \varrho } |$  denotes the volume of the \$d - d i m e n s i o n a l s p h e r e o f r a d i u s  $\varrho.$  F o r o u r o n e - d i m e n s i o n a l d o m a i n  $[0,\,L]$  ( i . e , \$d=1) we thus choose  $k _ { 2 } \, = \, 2 \, \varrho \, \lambda$  for the reaction-rate constant of the compartment-based model with sufficiently large compartment size  $h .$  Note that, on the other hand,  $h$  has to be small,  $h \ \ll \ L ,$  to ensure an appropriate level of spatial resolution (Erban and Chapman 2020).

<span id="page-16-0"></span>

## 3.2 Model Comparison

In contrast to the first-order reaction-diffusion system of Sect. 2, analytical insights for the dimerization system are relatively limited, so our studies solely rely on computational simulations. Figure 6 shows the steady-state gradients of molecules of A and  $B ,$  estimated from long-term stochastic simulations (using the Gillespie stochastic

sed -based -base

<span id="page-17-0"></span>Fig. 6 Steady-state solution for the morphogen system with the dimerization reaction. Comparison between Brownian dynamics (solid lines) and generalized mgRDME model (dots/ bars) for dimerization system (2). The parameters are given in (42), together with  $L =$  1,  $K _ { A } = 1 5$  and  $K _ { B } = 6 0$  (Color figure online)

simulation algorithm for the generalized mgRDME and temporal discretization for the Brownian dynamics model). The parameters are chosen as

<span id="page-17-1"></span> $k _ { 1 } = 5 0 , \ \ \lambda = 5 , \ \ \varrho = 0 . 0 2 , \ \ k _ { 2 } = 2 \, \lambda \, \varrho = 0 . 2 , \ \ k _ { 3 } = 2 , \ \ D _ { A } = 0 . 1 6 , \ \ D _ { B } = 0 . 0 1 .$  (42)

We observe a close agreement of the generalized compartment-based model with the particle-based dynamics for  $K _ { A } \, = \, 1 5$  and  $K _ { B } \, = \, 6 0 .$  In contrast, for the standard RDME the agreement is worse, as we will see in the next section.

## 3.2.1 Multi-Objective Optimization

For the spatial error in the distribution of molecules of B, we define the function  $\mathrm { e r r } ( K _ { B } , \, \gamma )$  in analogy to (36), that is, we choose the  $L ^ { 1 } - d i s t a n c e b e t w e e n t h e s t e a d y -$  state solutions of the compartment-based model and of the particle-based model. Let  $\bar { A } _ { \mathrm { t o t a l } }$  and  $\bar { B } _ { \mathrm { t o t a l } }$  be the total number of molecules of  $A$  and  $B ,$  respectively, of the compartment-based dimerization system at steady state. In contrast to the first-order reaction system of Sect. 2, we lack analytical expression for these quantities and can only estimate them via long-term numerical simulation. The cost function is defined similarly to (40), using the total propensity for diffusion and reaction events at steady state:

<span id="page-17-2"></span>which is quadratic in the total number of molecules of A at steady state. Figure 7 illustrates the error and cost values across different combinations of  $K _ { B }$  and  $\gamma .$  We see that the error is non-monotone in y, which differs from the behavior seen in Fig.5 for the first-order reaction-diffusion system. It seems that there is an optimal  $K _ { A } \approx 1 5 \, = \, 1 2 0 / 8 \, = \, 6 0 / 4 ,$  where the error in the spatial distribution of molecules of  $B$  is minimal. Next, we study the relationship between the parameter  $K _ { A }$  and the error in the A-population.

## Springer

|     |        |        |       |     | + standard                   |        |     |
|-----|--------|--------|-------|-----|------------------------------|--------|-----|
|     | 10 4   |        |       |     | generalized with K_{g}= 120  |        |     |
|     |        |        |       |     | O generalized with K_{0} =60 |        |     |
|     |        |        |       | *01 |                              |        |     |
|     | 160 3) |        |       |     |                              |        |     |
|     | 10 41  |        |       |     |                              |        |     |
|     |        | To0 5) |       |     |                              | “20 11 |     |
| 10  |        |        |       |     |                              |        |     |
| 0.2 | 0.3    | 0.4    | 0.5   | 0.6 | 0.7                          | 0.8    | 0.9 |
|     |        |        | error |     |                              |        |     |

<span id="page-18-0"></span> $\textbf { F i g . 7 }$  Error and cost for dimerization. Cost  $c ( K _ { B } , \, \gamma )$  as defined in (43) and error  $\mathrm { e r r } ( K _ { B } , \, \gamma )$  given by (36) for different values of  $K _ { B }$  and  $\gamma = K _ { B } / K _ { A }$  (given by labels of the form  $\{ K _ { B } \ \gamma$  ). The parameters are given in (42) and  $L =$  1. The results are calculated as averages over long-term simulations of length  $T = 5 0 0 0$  (Color figure online)

## 3.2.2 Error Analysis for Species A

Given the steady-state solution  $\bar { \bf A } = ( \bar { A } _ { 1 } , \, \bar { A } _ { 2 } , \, \ldots , \, \bar { A } _ { K _ { A } } )$  of the generalized mgRDME model for the dimerization system, we define

$$\bar { a } _ { h _ { A } } ( x ) : = \bar { A } _ { i } / h _ { A } \quad \mathrm { f o r } \, x \in [ ( i - 1 ) h _ { A } , i h _ { A } ) , \quad i = 1 , 2 , \ldots , K _ { A } , \qquad ( 4 4 )$$

as an approximation of the position-dependent steady-state concentration. Let further a(x) be the steady-state solution of the particle-based Brownian-dynamics dimerization model (these functions can be approximately determined via numerical simulation). We define the spatial error in A in analogy to (36):

<span id="page-18-2"></span><span id="page-18-1"></span>
$$\mathrm { ~ e r r _ { A } ( K _ { A } ) : = | \bar { a } - \tilde { a } _ { h } | | _ { L } | = \int _ { 0 } ^ { L } | \, \tilde { a } ( x ) \, | \, d x . \qquad \qquad ( 4 5 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

Moreover, we consider the difference in the total number of molecules of A at steady state:

As opposed to the first-order system, this quantity can here be non-zero. Figure 8 shows a non-monotone relationship between the number  $K _ { A }$  of compartments and the error erriotn1}(KA) in the total number of molecules of A. The non-monotony is due to a too small grid size  $h _ { A }$  for large  $\, K _ { A }$  which contrasts with the condition  $\varrho / h \ll 1$  necessary for a good approximation of second-order reactions, see also the paragraph on the relation between the rates on page 17. The spatial error err4(K4) in  $A ,$  however, decreases monotonically with increasing  $K _ { A }$  because a higher level of spatial resolution dominates the effect of losing second-order reactions for the parameter values under consideration. The error in  $B$  (also depicted in Fig. 8b) results from a combination of the two errors in  $A$  and minimizes for  $K _ { A } \approx$  16, which is in consistency with the observations from Fig. 7. The standard compartment-based model would mean to choose  $K _ { A } = K _ { B } = 6 0 ,$  inducing a comparatively large error both in

| + L | al |
|-----|----|
|-----|----|

 $K _ { A }$  25

<span id="page-19-1"></span>(a) Error in  $A$  for  $K _ { A } \in \{ 2 , 4 , . . . , 2 0 \}$ 

(b) Error in  $A$  and  $B$  for  $K _ { A } \in \{ 1 0 , 1 2 , . . . , 2 0 \}$ 

 $\textbf { F i g . 8 }$  Error in dependence on  $K _ { A } \mathbf { s }$  Spatial error  $\mathrm { e r r } _ { A } ( K _ { A } )$  given by (45) (blue) and difference  $\mathrm { e r r _ { t o t a l } ( K _ { A } ) }$  in the total number of molecules of  $A$  defined in (46) (red) for different values of  $K _ { A } . \mathbf { b }$  zoom-in, additionally contains the spatial error in  $B , \exp ( K _ { B } , \, \gamma )$  defined in (36), for fixed  $K _ { B } = 6 0$  and  $\gamma = K _ { B } / K _ { A } .$  Parameters given in (42),  $L = 1 .$  Long-term simulation of length  $T = 1 0 ^ { 5 }$  (Color figure online)

## the steady-state distribution of molecules of $B$ and in the total number of molecules of A.

We conclude that the generalized mgRDME model clearly outperforms the standard RDME model for the fast-slow morphogen gradient system with dimerization. For practical implementation, the preliminary choice of  $\gamma = \sqrt { D _ { A } / D _ { B } }$  must be carefully adapted to ensure the adequate reproduction of the second-order reaction  $A + A \to B .$  Specifically, the compartment size  $h _ { A }$  should be chosen based on the reaction radius  $\varrho$  and the rate constant \lambda of the particle-based Brownian dynamics model, in conjunction with the diffusion constant  $D _ { A } ,$  but independently of the spatial resolution for  $B .$  This stands in clear contrast to the first-order chemical reaction system discussed in Sect. 2, where it is more appropriate to adjust the compartment sizes of the two species relative to each other.

Remark 3 (Second-order moments for dimerization) Unlike the first-order reactiondiffusion model discussed in Sect. 2, the variance of the particle numbers does here not align with the mean. Consequently, one could also take the variance as a value for comparing the models and defining an error. A numerical comparison of mean and variance for the different model types in application to the dimerization system is given in Appendix 5.3.

<span id="page-19-0"></span>

## 4 Discussion

We examined the mgRDME as a generalized compartment-based stochastic model for reaction-diffusion kinetics, allowing adapted compartment sizes for different chemical species involved in the system (Cao and Erban 2014). When employed in the context of morphogen gradient formation, the mgRDME is useful because the signaling molecules A typically diffuse at a significantly higher rate than the morphogen molecules B, and this diffusion rate is a deciding factor for the choice of compartment sizes. In comparison to the ‘ground truth’particle-based Brownian dynamics simulations, the mgRDME outperforms the standard RDME in terms of accuracy and

reduction in numerical cost, both for the first-order reaction network (1) and for the considered dimerization system (2).

Interestingly, however, the structural dependence of error and cost on the compartment sizes is different for the two settings. In Sect. 2, the optimal ratio between grid sizes hinges on the relationship between the diffusion rates. Conversely, in the system involving second-order dimerization in Sect. 3, an optimal grid size for species A can be identified (in dependence on its diffusion constant) regardless of the grid size allocated to  $B .$  In particular, the approximation error may be reduced for the dimerization system by choosing a smaller compartment size for  $B$  while maintaining the fixed optimal compartment size for A undergoing dimerization. This goes in line with the observations for standard RDMEs, where very small grid sizes result in the elimination of the second-order reactions in three-dimensional simulations. The mgRDME circumvents this issue by allowing larger compartment sizes for species involved in second-order reactions. At the same time, the increase in numerical cost caused by the higher spatial resolution is limited, since it is selectively applied solely to the slow species.

Our work shows that the mgRDME is a valuable framework for studying reactiondiffusion systems with multiple scales in diffusion speed of involved molecules. The multi-grid approach enables the mgRDME to more effectively accommodate the inherent conditions of the system under consideration, surpassing the adaptability of the standard RDME. The models with different diffusion constants naturally take into account variations in the sizes of simulated biomolecules, which can range from small calcium ions to relatively large vesicles in neurotransmission dynamics (Ernst et al. 2022, 2023), or from relatively small  $G - a c t i n$  monomers to larger  $F - a c t i n f i l a m e n t s i n$  simulations of actin dynamics (Zhuravlev et al. 2010; Dobramysl et al. 2016b). Other applications of mgRDME include pattern formation based on Turing instability which requires different diffusion constants (Turing 1952).

In this work, we have solely considered dynamics in one-dimensional domains – a restriction which is well justified in the context of applications to developmental biology, where a modeller is interested in the pattern formation along one body axis, see, for example, dorsal-ventral patterning studies of the Drosophila embryo (Shimmi and O’Connor 2003; Shimmi et al. 2005). However, the ideas behind mgRDME can be extended to higher-dimensional models, where unstructured meshes are used to discretize complex geometries of biological problems to better resolve curved inner and outer boundaries than by Cartesian meshes (Engblom et al. 2009; Isaacson and Zhang 2018). The mgRDME framework can use different unstructured meshes to describe different chemical species. In two-dimensional domains, voxels are typically triangles or squares, while in three-dimensional domains, they are tetrahedra or cubes (Engblom et al. 2009; Earnest et al. 2015). The reaction rates in mgRDME are constructed based on the shared area or volume between voxels from different grids, as described in Remark 1, where we note that voxels may partially overlap, see Eq. (25). The relative size of different meshes (parameterized by γ in this paper) can be chosen to optimize the number of spatial diffusive jumps for species with varying diffusion constants. Additionally, different species can be modeled at different resolutions to improve simulation accuracy, with voxel sizes optimized to enhance the description of bimolecular reactions. In applications where high accuracy is only needed in a subset

of the simulation domain (Dobramysl et al. 2016a, b), unstructured meshes also offer the advantage to refine the mesh in the corresponding subset of the computational domain. In particular, to model complex biological systems, the mgRDME approach (of using different meshes for different chemical species) can be combined with using different mesh resolutions (for the same chemical species) in different parts of the computational domain (Kang and Erban 2019).

While our analysis has been based on comparison of steady-state morphogen gradients, another way to compare the RDME and mgRDME could be by comparing the first collision and mean reaction times in these models (Li et al. 2018), which can be achieved by analyzing the random walks on lattices (Montroll 1969). Also a comparison of time-dependent solutions (in contrast to steady-state distributions investigated here) is another interesting topic for future research on mgRDME models of morphogen gradients (Saunders and Howard 2009). Some approaches to compare the dynamics prior to steady state are discussed in Appendix 5.4.

<span id="page-21-0"></span>

## 5 Appendix

## 5.1 Time-Dependent Solution of the Brownian Dynamics Model

Let us consider that initially there are no molecules of  $A$  and  $B$  in the system, i.e.

<span id="page-21-1"></span>
$$a ( x , 0 ) \equiv b ( x , 0 ) \equiv 0 . \qquad \qquad \qquad \qquad ( 4 7 )$$

To solve the PDE (27) with initial condition (47) and boundary conditions (29), we apply the Fourier transform

$$\hat { a } ( \xi , t ) = \int _ { - \infty } ^ { \infty } a ( x , t ) \, e ^ { - i \xi x } \, \mathrm { d } x .$$

We get

∂ a ∂t (ξ , t) =−(DA ξ 2 + k2) a(ξ , t) + 2k1.

Solving this ODE with the initial condition  $\widehat { a } ( \xi ,$  0)  $= 0 ,$  we obtain

<span id="page-21-2"></span>a(ξ , t) = 2k1 DA ξ 2 + k2 − 2k1 exp DA ξ 2 + k2 − - DA ξ 2 + k2  t ! . (48)

Since the Fourier transform of the convolution is the product of Fourier transforms, we can calculate the Fourier inverse of (48) as follows

a(x, t) =¯a(x) − exp[−k2t] 2 √DAtπ −∞ ∞ a¯(y) exp" −(x − y)2 4DAt # dy,

where a¯(x) is given by (33). In total we get

## a(x, t) = k1 √DA k2 exp − DA k2 |x|

− exp[−k2t] 2 √DAtπ −∞ ∞ √DA k2 k1 exp − DA k2 |y| exp" −(x − y)2 4DAt # dy

<span id="page-22-2"></span>(49)

for the average number of A-molecules in the interval [x, x + dx] at time t.

<span id="page-22-1"></span>

## 5.2 ’Ground Truth’ Solution: The Effect of the Boundary

Here, we derive the ‘ground truth’ solution of Sect. 2.3 for zero-flux boundary conditions. To do that, we will solve the ODEs (30–31) on the interval [−L, L] with zero-flux boundary conditions

$$\frac { \mathrm { d } \bar { x } } { \mathrm { d } x } ( - L ) = \frac { \mathrm { d } \bar { a } } { \mathrm { d } x } ( L ) = \frac { \mathrm { d } \tilde { b } } { \mathrm { d } x } ( L ) = 0 . \qquad \qquad ( 5 0 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

Then the solution to the steady-state Eq. (30) is

a¯(x) = k1 √DA k2 exp − DA k2 |x| 1 + exp\$ 1 − exp[−2L 2 % k2 DA √k2/DA] ,  |x|− L &

## which for L →∞ agrees with (33), independently of the parameter values.

<span id="page-22-3"></span>

## 5.3 Second-Order Moments for Dimerization

Means and variances of the numbers of molecules of  $A$  and  $B$  for the dimerization system are plotted in Fig. 9, sampled from long-term simulations of the (standard) compartment-based model and of the particle-based Brownian dynamics model, choosing boxes of size  $h \, = \, 1 / 1 5 .$  We observe a good agreement of the two model types.

<span id="page-22-0"></span>

## 5.4 Dynamics Prior to Steady State

To study the dynamics prior to steady state we consider the average positions of the particles in dependence on time. These investigations are done for the first-order morphogen gradient system introduced in Sect. 2.

51

• A sample mean compartment-based ]A sample mean particle-based

• B sample mean compartment-based  $\mathbb { B }$  sample mean particle-based

0.4 0.6 position 0.8 (a) mean A  $. 2$   $0 . 4$  0.6 position 0.8 (b) mean  $B$ 

A sample variance compartment-based ]A sample variance particle-based

• B sample variance compartment-based  $\exists B$  sample variance particle-based

0.6 position 0.8

(c) variance  $A$ 

<span id="page-23-3"></span><span id="page-23-2"></span>osition 0.8 (d) variance  $B$ 

<span id="page-23-0"></span>Fig. 9 Comparison of mean and variance for dimerization. Mean/variance of the number of molecules of  $\textbf { a / c }$  species \$A \ \$ and  $\mathbf { b / d }$  species \$B \ \$ at steady state estimated by long-term simulation \$( T = 2 0 0 0 ) \ \$ of the compartment-based model (red circles) and of the particle-based Brownian dynamics model (grey bars). Rate constants given in (42),  $L =$  1,  $K _ { A } = K _ { B } = 1 5 ~ ( \gamma = 1 , h = 1 / 1 5 )$  (Color figure online)

## 5.4.1 Average Positions for ‘Ground Truth’ Solution

## We define the moments

$$\langle 1 \rangle _ { a } ( t ) : = \int _ { 0 } ^ { \infty } a ( x , t ) \, \mathrm { d } x \, , \qquad \langle x \rangle _ { a } ( t ) : = \int _ { 0 } ^ { \infty } x \, a ( x , t ) \, \mathrm { d } x \, , \qquad ( 5 1 )
\int _ { 0 } ^ { \infty }$$

Multiplying Eqs. (27-28) by  $1$  and  $x ,$  respectively, and integrating over  $x \, \in [ 0 , \, \infty ) ,$  we obtain the ODEs

<span id="page-23-1"></span> $\frac { \mathrm { d } ( 1 ) _ { a } } { \mathrm { d } t } = - k _ { 2 } \, ( 1 ) _ { a } \, + \, k _ { 1 } ,$  (53)

= D, a (0.1) - 2 (x).

(54)

$$\frac { \mathrm { d } ( 1 ) _ { b } } { 4 . } = k _ { 2 } \left ( 1 \right ) _ { a } \, - \, k \, \, \, \qquad \qquad \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

\$100 ={0}(0.0)+ 2 (x) = (x).

(56)

<span id="page-24-0"></span>Substituting  $x = 0$  into Eq. (49), we have

(57)

Substituting (57) into Eq. (54), we get

Using the initial condition  $\langle x \rangle _ { a } ( 0 ) = 0 ,$  we obtain

<span id="page-24-3"></span>
$$\ \langle x \rangle _ { a } ( t ) = \exp [ - k _ { 2 } \, t ] \int _ { 0 } ^ { t } \exp \left [ k _ { 2 } \tau \right ] \frac { k _ { 1 } \sqrt { b _ { A } } } { \sqrt { k _ { 2 } \tau } } \right ] \mathrm { d } \tau \ .$$

<span id="page-24-1"></span>Using the initial condition  $\langle 1 \rangle _ { a } ( 0 ) = 0$  in Eq. (53), we obtain

<span id="page-24-2"></span>for the average total number of A-particles at time t. Consequently, the average position  $x _ { a } ( t )$  of A-particles as a function of time  $t > 0$  is

We note that this ratio is independent of our choice of the rate constant k1.

Substituting (59) into Eq. (55), we get

2

x n(t)  $x _ { b } ( t )$ 

time

(a) Average number of particles

2 3 time 4 5

(b) Average positions

<span id="page-25-2"></span>Fig. 10 ‘Ground truth’ dynamics prior to steady state. a Average total number of particles (1)}(t) and  $\langle 1 \rangle _ { b } ( t )$  given by (59) and (61), respectively.  $\mathbf { b }$  Average positions  $x _ { a } ( t )$  and  $x _ { b } ( t )$  depending on time given by (60) and (62). Rate constants given in (18) (Color figure online)

Using the initial condition  $\langle 1 \rangle _ { b } ( 0 ) = 0 ,$  we obtain

<span id="page-25-0"></span>

## (1) st) = epl-ks.) [ ep[ks r]k: (1 – epl-kr) 4r

$$= \frac { k _ { 1 } } { k _ { 3 } } + \frac { k _ { 1 } k - k _ { 3 } } { ( 4 - k _ { 3 } ) k _ { 3 } } + \frac { k _ { 1 } \cos k _ { 2 } } { k _ { 3 } - k _ { 3 } } \qquad \qquad ( 6 1 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

for the average total number of \$B - p a r t i c l e s a t t i m e  $t.$  S u b s t i t u t i n g ( 5 8 ) i n t o E q . ( 5 6 ) ,\$ we get

<span id="page-25-1"></span>The average position  $x _ { b } ( t )$  of \$B - p a r t i c l e s a t t i m e  $t$  i s g i v e n b y\$

Fig. 10 shows the average total number of particles and their average positions depending on time for the ‘ground truth’ model.

## 5.4.2 Time-Dependent Solution for the Compartment-Based Model

In analogy to the densities  $a ( x , \, t )$  and  $b ( x , \, t )$  defined for the particle-based model, let  $a _ { i } ( t )$  and  $b _ { i } ( t )$  denote the average number of particles of species  $A$  and  $B$  in box  $i$  at time  $t ,$  respectively. We set  $\mathbf { a } ( t ) = ( a _ { i } ( t ) ) _ { i = 1 , \ldots , K _ { A } }$  and  $\mathbf { b } ( t ) = \, ( b _ { i } ( t ) ) _ { i = 1 , \ldots , K _ { B } }$  and

Springer

| 0.07    |                                         | 0.07   |                                                                                                          |
|---------|-----------------------------------------|--------|----------------------------------------------------------------------------------------------------------|
|         | x (t)                                   |        | \$- x _ { B } ^ { ( 1 ) } ( t ) - x _ { b } ( t ) , \, K _ { B } = 1 2 0 \big  \$                        |
| § 0.06  | \$x _ { b } ( t )\$                     | \$0.06 | \$- x _ { B } ^ { ( 1 ) } ( t ) - x _ { b } ( t ) , \ K _ { B } = 1 0 0 \big  \$                         |
|         | ax{0}(t)                                |        | \$- x _ { B } ^ { ( 1 ) } ( t ) - x _ { b } ( t ) , \, K _ { B } = 8 0 \, \,\$                           |
| 2. 0.05 |                                         | 0.05   |                                                                                                          |
|         | \$- x _ { b } ( t ) \big  \$            |        | \$- x _ { B } ^ { ( 1 ) } ( t ) - x _ { b } ( t ) , \, K _ { B } = 6 0\$                                 |
| = 0.04  | \$- \, x _ { b } ( t ) \big  \$ x{0}(t) | 0.04   | \$- x _ { B } ^ { ( 1 ) } ( t ) - x _ { b } ( t ) , \ K _ { B } = 4 0 \ \$ ax{0}{3}(t)-x_{0}(t), K{0}=20 |
| 0.03    |                                         | - 0.03 |                                                                                                          |
| 0.0     |                                         | § 0.02 |                                                                                                          |
|         |                                         |        | 4 3 5                                                                                                    |
|         | time                                    |        | \$\mathrm { \ t i m e }\$                                                                                |

(a) Difference depending on  $\gamma$  (b) Difference depending on  $K _ { B }$   $\ ( K _ { B } = 1 2 0 )$   $\ ( \gamma = 1 )$ 

<span id="page-26-4"></span>(t) Fig. 11 Time-dependent differences of average positions. a Differences between average positions \$x of mgRDME, see (69), and average positions x0(t) of the ‘ground truth’ solution, see (62), for various y and fixed  $K _ { B } = 1 2 0 . \, \mathbf { b }$  Differences between average positions of the RDME (fixed  $\gamma = 1 )$  and the  $^ { * }$  ground truth’ solution for various  $K _ { A } = K _ { B } .$  Rate constants given in (18) (Color figure online)

have

<span id="page-26-1"></span><span id="page-26-0"></span>1 5-24. )

$${ \frac { \mathrm { d } b } { \mathrm { d } t } } = \left ( { \frac { D _ { p } } { h _ { \gamma } } } S _ { B } - k _ { 3 } I _ { B } \right ) { \mathrm { ~ d ~ } } } , \qquad \qquad \qquad ( 6 4 ) \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \; \;$$

for  $I _ { A } , S _ { A } \in \mathbb { R } ^ { K _ { A } \times K _ { A } } , \ I _ { B } , S _ { B } \ \in \mathbb { R } ^ { K _ { B } \times K _ { B } } , \ \mathbf { e } _ { 1 } }$  and  $M \in \mathbb { R } ^ { K _ { A } }$  defined in Sect. 2.2. We note that the standard RDME is included by setting  $K _ { A } = K _ { B } .$  The ODEs (63)-(64) are the analogue to Eqs. (27–28) of the ‘ground truth’ dynamics. Together with  $\mathbf { a } ( 0 ) = 0$  we obtain

<span id="page-26-2"></span>
$$\mathbf { a } ( t ) = W _ { A } ^ { - 1 } \left [ \exp { ( W _ { A } t ) } - I _ { A } \right ] k _ { 1 } \qquad \qquad \qquad \qquad ( 6 5 )$$

for

<span id="page-26-3"></span>
$$\ W _ { A } : = \frac { D _ { A } } { \hbar A } S _ { A } - k _ { 2 } I _ { A } \in \mathbb { R } ^ { K _ { A } \times K _ { A } } .$$

Moreover, using  $\mathbf { b } ( 0 ) = 0$  and defining  $W _ { B } : = ( D _ { B } / h _ { B } ^ { 2 } ) S _ { B } - k _ { 3 } I _ { B } \in \mathbb { R } ^ { K _ { B } } ,$  we obtain (Perko 2013)

for the vector of time-dependent average numbers of  $B$  -particles in the boxes at time  $t .$  We have limt\infty\_{0}\$ a(t) = A and limt,{\infty} b(t) = B for the steady-state solutions A, B given in (16) and (24).

## 5.4.3 Average Positions for the Compartment-Based Model

<span id="page-27-2"></span><span id="page-27-1"></span>In analogy to the moments defined in (51)-(52), we set

In the definitions of  $\langle x \rangle _ { A } ( t )$  and  $\langle x \rangle _ { B } ( t ) ,$  we choose the centers  $( i \mathrm { ~ - ~ } 0 . 5 ) h _ { A }$  and  $( i \mathrm { ~ - ~ } 0 . 5 ) h _ { B }$  of the boxes for the spatial value, which is justified by the assumption that the particles in each box are spatially well mixed so that their positions follow a uniform distribution.

By summing Eqs. (63-64) over the discrete index  $i$  we obtain the same ODEs for the moments (1) 4(t) and (1) (t) as in the ‘ground truth’model; see Eqs. (53) and (55). This means that the average total number of particles  $\langle 1 \rangle _ { A } ( t )$  and  $\langle 1 \rangle _ { B } ( t )$  is again given by Eqs. (59) and (61). We calculate the moments  $\langle x \rangle _ { A } ( t )$  and  $\langle x \rangle _ { B } ( t )$  according to their definition in (67)-(68), using the solutions  $\mathbf { a } ( t )$  and  $\mathbf { b } ( t )$  given in (65) and (66). Given  $\gamma = K _ { B } / K _ { A } ,$  the mean positions  $x _ { A } ^ { ( \gamma ) } ( t )$  and  $x _ { B } ^ { ( \nu ) } ( t )$  of  $A$  and  $B$  molecules at time  $t$  are

<span id="page-27-0"></span>

## in analogy to (60) and (62), respectively.

Figure 11 shows the time-dependent differences between the average positions  $x _ { B } ^ { \mathrm { ~ V ~ } } ( t )$  of the mgRDME and the average positions  $x _ { b } ( t )$  of the  $^ { \cdot }$  ground truth’ solution for various  $\gamma$  and  $K _ { B } .$  We observe that for every time  $t >$  0, the error increases with the ratio  $\gamma = K _ { B } / K _ { A }$  for fixed  $K _ { B } = 1 2 0 ,$  see Fig. 11a. Similarly, the error increases with decreasing values of  $K _ { B } ,$  see Fig. 11b. In contrast to the results from Sect.2.4, the approximation quality is here more robust against changes in  $K _ { B } ;$  this can be explained by the fact that the reference value is given by the average positions and not by the spatial distribution of molecules of B.

The mean position  $x _ { B } ^ { \mathrm { ~ V ~ } } ( t )$  of the compartment-based model always exceeds the mean position \$x{0}(t) of the ‘ground truth’ solution. Especially for the initial time phase (small t), the difference is large. This is due to the fact that the production of molecules of  $A$  is located at position  $x \, = \, 0 .$  In the mgRDME, the newly produced particles immediately get well mixed in the first box I0,  $h _ { A } ] ,$  while in the  $^ { - }$  ground truth’model it takes some time for the particles to spread out in space. The mixing in the mgRDME implies a shift of the average to the right, which has a larger effect especially in the beginning, where the population is small and mostly located near the

origin. This effect, which also propagates to  $B$  via conversion, can be interpreted as an artificial dispersion induced by the mgRDME.

Author Contributions Both authors formulated the problem, carried out the theoretical analysis of the models and wrote the paper. SW carried out illustrative numerical simulations.

Funding This work was supported by the Engineering and Physical Sciences Research Council, grant number EP/V047469/1, awarded to Radek Erban. We further acknowledge the support of Deutsche Forschungsgemeinschaft through CRC 1114 and through Germany’s Excellence Strategy–The Berlin Mathematics Research Center MATH+ (EXC-2046/1, project 390685689), awarded to Stefanie Winkelmann.

Data Availability In compliance with EPSRC’s open access initiative, the data in this paper is available from: https://doi.org/10.5281/zenodo.11072450.

## Declarations

Conflict of interest Authors have no conflict of interest.

Ethical Approval Not applicable.

Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and indicate if changes were made. The images or other third party material in this article are included in the article’s Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is not included in the article’s Creative Commons licence and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this licence, visit http://creativecommons.org/licenses/by/4.0/.

<span id="page-28-10"></span><span id="page-28-9"></span><span id="page-28-8"></span><span id="page-28-7"></span><span id="page-28-6"></span><span id="page-28-5"></span><span id="page-28-4"></span><span id="page-28-3"></span><span id="page-28-2"></span><span id="page-28-1"></span><span id="page-28-0"></span>

## References

Agbanusi IC, Isaacson SA (2014) A comparison of bimolecular reaction models for stochastic reactiondiffusion systems. Bull Math Biol 76(4):922–946 Andrews SS (2012) Spatial and stochastic cellular modeling with the Smoldyn simulator. Methods and Protocols, Bacterial Molecular Networks, pp 519–542 Andrews SS, Addy NJ, Brent R, Arkin AP (2010) Detailed simulations of cell biology with Smoldyn 2.1. PLoS Comput Biol 6(3):1000705 Ben-Avraham D, Redner S (1986) Kinetics of  $n - s p e c i e s a n n i h i l a t i o n : m e a n - f i e l d a n d d i f f u s i o n - c o n t r o l l e d$  limits. Phys Rev A 34(1):501 Bergmann S, Sandler O, Sberro H, Shnider S, Schejter E, Shilo B, Barkai N (2007) Pre-steady-state decoding of the Bicoid morphogen gradient. PLoS Biol 5(2):46 Cao Y, Erban R (2014) Stochastic Turing patterns: analysis of compartment-based approaches. Bull Math Biol 76:3051–3069 Chatterjee A, Vlachos D, Katsoulakis M (2004) Spatially adaptive lattice coarse-grained Monte Carlo simulations for diffusion of interacting molecules. J Chem Phys 121:11420–11431 Dai J, Seider W, Sinno T (2008) Coarse-grained lattice kinetic Monte Carlo simulation of systems of strongly interacting particles. J Chem Phys 128:194705 Denizot A, Arizono M, Nägerl U, Soula H, Berry H (2019) Simulation of calcium signaling in fine astrocytic processes: effect of spatial properties on spontaneous activity. PLoS Comput Biol 15(8):1006795 Dobramysl U, Rüdiger S, Erban R (2016) Particle-based multiscale modeling of calcium puff dynamics. Multiscale Modelling Simulation 14(3):997–1016 Dobramysl U, Papoian G, Erban R (2016) Steric effects induce geometric remodeling of actin bundles in filopodia. Biophys J 110:2066–2075

Doi M (1976) Stochastic theory of diffusion-controlled reaction. J Phys A: Math Gen 9(9):1479 Doi M (1976) Second quantization representation for classical many-particle system. J Phys A: Math Gen 9(9):1465 Earnest TM, Lai J, Chen K, Hallock MJ, Williamson JR, Luthey-Schulten Z (2015) Toward a whole-cell model of ribosome biogenesis: kinetic modeling of SSU assembly. Biophys J 109(6):1117–1135 Engblom S, Ferm L, Hellander A, Lötstedt P (2009) Simulation of stochastic reaction-diffusion processes on unstructured meshes. SIAM J Sci Comput 31(3):1774–1797 Erban R, Chapman SJ (2007) Reactive boundary conditions for stochastic simulations of reaction-diffusion processes. Phys Biol 4(1):16–28 Erban R, Chapman SJ (2009) Stochastic modelling of reaction-diffusion processes: algorithms for bimolecular reactions. Phys Biol 6(4):046001 Erban R, Chapman SJ (2020) Stochastic Modelling of Reaction-Diffusion Processes. Cambridge University Press, Cambridge Erban R, Flegg M, Papoian G (2014) Multiscale stochastic reaction-diffusion modelling: application to actin dynamics in filopodia. Bull Math Biol 76(4):799–818 Ernst A, Schütte C, Sigrist SJ, Winkelmann S (2022) Variance of filtered signals: characterization for linear reaction networks and application to neurotransmission dynamics. Math Biosci 343:108760 Ernst A, Unger N, Schütte C, Walter AM, Winkelmann S (2023) Rate-limiting recovery processes in neurotransmission under sustained stimulation. Math Biosci 362:109023 Fange D, Elf J (2006) Noise-induced Min phenotypes in E. coli. PLoS Comput Biol 2(6):637–648 Flegg M, Chapman J, Erban R (2012) The two-regime method for optimizing stochastic reaction-diffusion simulations. J R Soc Interface 9(70):859–868 Flegg M, Hellander S, Erban R (2015) Convergence of methods for coupling of microscopic and mesoscopic reaction-diffusion simulations. J Comput Phys 289:1–17 Gillespie DT (1977) Exact stochastic simulation of coupled chemical reactions. J Phys Chem 81(25):2340– 2361 Gillespie DT (2007) Stochastic simulation of chemical kinetics. Annu Rev Phys Chem 58:35–55 Gregor T, Wieschaus E, McGregor A, Bialek W, Tank D (2007) Stability and nuclear dynamics of the bicoid morphogen gradient. Cell 130(1):141–152 Hellander S, Petzold L (2016) Reaction rates for a generalized reaction-diffusion master equation. Phys Rev E 93(1):013307 Hellander S, Hellander A, Petzold L (2012) Reaction-diffusion master equation in the microscopic limit. Phys Rev E 85(4):042901 Hellander S, Hellander A (2020) Hierarchical algorithm for the reaction-diffusion master equation. Journal of Chemical Physics 152(3) Isaacson SA (2013) A convergent reaction-diffusion master equation. Journal of Chemical Physics 139(5) Isaacson SA (2008) Relationship between the reaction-diffusion master equation and particle tracking models. J Phys A: Math Theor 41(6):065003 Isaacson SA (2009) The reaction-diffusion master equation as an asymptotic approximation of diffusion to a small target. SIAM J Appl Math 70(1):77–111 Isaacson SA, Peskin CS (2006) Incorporating diffusion in complex geometries into stochastic chemical kinetics simulations. SIAM J Sci Comput 28(1):47–74 Isaacson SA, Zhang Y (2018) An unstructured mesh convergent reaction-diffusion master equation for reversible reactions. J Comput Phys 374:954–983 Iyer K, Prabhakara C, Mayor S, Rao M (2023) Cellular compartmentalisation and receptor promiscuity as a strategy for accurate and robust inference of position during morphogenesis. Elife 12:79257 Jahnke T, Huisinga W (2007) Solving the chemical master equation for monomolecular reaction systems analytically. J Math Biol 54:1–26 Kang H, Erban R (2019) Multiscale stochastic reaction-diffusion algorithms combining Markov chain models with stochastic partial differential equations. Bull Math Biol 81(8):3185–3213 Kang K, Redner S (1984) Scaling approach for the kinetics of recombination processes. Phys Rev Lett 52(12):955 Kang H, Zheng L, Othmer H (2012) A new method for choosing the computational cell in stochastic reaction-diffusion systems. J Math Biol 65:1017–1099 Kicheva A, Pantazis P, Bollenbach T, Kalaidzidis Y, Bittig T, Julicher F, González-Gaitán M (2007) Kinetics of morphogen gradient formation. Science 315(5811):521–525

Kolomeisky AB (2011) Formation of a morphogen gradient: acceleration by degradation. J Phys Chem Lett 2(13):1502–1505 Li F, Chen M, Erban R, Cao Y (2018) Reaction time for trimolecular reactions in compartment-based reaction-diffusion models. J Chem Phys 148(20):204108 Li F, Cao Y (2012) Multiscale discretization for reaction diffusion systems. In: Proceedings of the International Conference on Bioinformatics & Computational Biology (BIOCOMP), p. 242 Lipkova J, Zygalakis K, Chapman J, Erban R (2011) Analysis of Brownian dynamics simulations of reversible bimolecular reactions. SIAM J Appl Math 71(3):714–730 Montefusco A, Schütte C, Winkelmann S (2023) A route to the hydrodynamic limit of a reaction-diffusion master equation using gradient structures. SIAM J Appl Math 83(2):837–861 Montroll E (1969) Random walks on lattices. III. calculation of first-passage times with application to exciton trapping on photosynthetic units. J Math Phys 10(4):753–765 Murray J (2002) Math Biol II: Spatial Models Biomed Appl. Springer, New York Perko L (2013) Differential Equations and Dynamical Systems, vol 7. Springer, New York Rogers K, Schier A (2011) Morphogen gradients: from generation to interpretation. Annu Rev Cell Dev Biol 27(1):377–407 Saunders T, Howard M (2009) When it pays to rush: interpreting morphogen gradients prior to steady-state. Phys Biol 6(4):046020 Shimmi O, O’Connor M (2003) Physical properties of Tld, Sog, Tsg and Dpp protein interactions are predicted to help create a sharp boundary in Bmp signals during dorsoventral patterning of the Drosophila embryo. Development 130(19):4673–4682 Shimmi O, Umulis D, Othmer H, Connor M (2005) Faciliated transport of a Dpp/Scw heterodimer by Sog/Tsg leads to robust patterning of the Drosophila blastoderm embryo. Cell 120(6):873–886 Teimouri H, Kolomeisky AB (2018) Discrete-state stochastic modeling of morphogen gradient formation. Methods and Protocols, Morphogen Gradients, pp 199–221 Turing A (1952) The chemical basis of morphogenesis. Philosophical Trans R Soc Lond 237:37–72 Winkelmann S, Schütte C (2020) Stochastic Dynamics in Computational Biology, vol 645. Springer, Cham Winkelmann S, Zonker J, Schütte C, Conrad ND (2021) Mathematical modeling of spatio-temporal population dynamics and application to epidemic spreading. Math Biosci 336:108619 Winkelmann S, Schütte C (2016) The spatiotemporal master equation: Approximation of reaction-diffusion dynamics via Markov state modeling. Journal of Chemical Physics 145(21) Wolpert L, Tickle C, Arias AM (2015) Principles of Development. Oxford University Press, Oxford Zhuravlev P, Der B, Papoian GA (2010) Design of active transport must be highly intricate: a possible role of myosin and ena/VASP for G-actin transport in filopodia. Biophys J 98:1439–1448

<span id="page-30-18"></span><span id="page-30-17"></span><span id="page-30-16"></span><span id="page-30-15"></span><span id="page-30-14"></span><span id="page-30-13"></span><span id="page-30-12"></span><span id="page-30-11"></span><span id="page-30-10"></span><span id="page-30-9"></span><span id="page-30-8"></span><span id="page-30-7"></span><span id="page-30-6"></span><span id="page-30-5"></span><span id="page-30-4"></span><span id="page-30-3"></span><span id="page-30-2"></span><span id="page-30-1"></span><span id="page-30-0"></span>Publisher’s Note Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.
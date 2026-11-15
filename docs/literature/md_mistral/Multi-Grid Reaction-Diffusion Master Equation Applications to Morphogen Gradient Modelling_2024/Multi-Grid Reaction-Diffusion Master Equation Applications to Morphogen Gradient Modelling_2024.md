Bulletin of Mathematical Biology (2025) 87:6 https://doi.org/10.1007/s11538-024-01377-y

[No text]

Society for Mathematical Biology

ORIGINAL ARTICLE

![img-0.jpeg](img-0.jpeg)

# Multi-Grid Reaction-Diffusion Master Equation: Applications to Morphogen Gradient Modelling

Radek Erban \( ^{1} \)  · Stefanie Winkelmann \( ^{2} \)

Received: 3 May 2024 / Accepted: 21 October 2024 / Published online: 27 November 2024  
© The Author(s) 2024

## Abstract

The multi-grid reaction-diffusion master equation (mgRDME) provides a generalization of stochastic compartment-based reaction-diffusion modelling described by the standard reaction-diffusion master equation (RDME). By enabling different resolutions on lattices for biochemical species with different diffusion constants, the mgRDME approach improves both accuracy and efficiency of compartment-based reaction-diffusion simulations. The mgRDME framework is examined through its application to morphogen gradient formation in stochastic reaction-diffusion scenarios, using both an analytically tractable first-order reaction network and a model with a second-order reaction. The results obtained by the mgRDME modelling are compared with the standard RDME model and with the (more detailed) particle-based Brownian dynamics simulations. The dependence of error and numerical cost on the compartment sizes is defined and investigated through a multi-objective optimization problem.

Keywords Reaction-diffusion master equation · Multi-grid methods · Morphogen gradient formation · Stochastic simulation algorithms

## 1 Introduction

Compartment-based stochastic reaction-diffusion models have been used for modelling a range of biological processes, including Min-protein oscillations (Fange and Elf 2006) and ribosome biogenesis (Earnest et al. 2015) in E. coli, calcium signalling (Denizot et al. 2019; Dobramysl et al. 2016a), gene expression (Winkelmann and Schütte 2016), actin transport in filopodia (Zhuravlev et al. 2010; Erban et al.

Radek Erban
erban@maths.ox.ac.uk

Stefanie Winkelmann
winkelmann@zib.de

\( ^{1} \)  Mathematical Institute, University of Oxford, Radcliffe Observatory Quarter, Woodstock Road, Oxford OX2 6GG, UK

\( ^{2} \)  Zuse Institute Berlin (ZIB), Takustrasse 7, 14195 Berlin, Germany

Springer

6 Page 2 of 31

R. Erban, S. Winkelmann

2014) and epidemic spreading (Winkelmann et al. 2021). They can be simulated using algorithms for continuous-time discrete-space Markov chains and are mathematically described using the reaction-diffusion master equation (RDME) for the probability mass function (Erban and Chapman 2020; Winkelmann and Schütte 2020). At the analytical level, the RDME builds a bridge between microscopic models for reaction-diffusion processes and macroscopic approximations through partial differential equations (Flegg et al. 2012; Montefusco et al. 2023; Flegg et al. 2015; Kang and Erban 2019).

To formulate a standard compartment-based model, the computational domain is discretized into compartments and diffusion is modelled as a jump process between the compartments (Winkelmann and Schütte 2020). Considering three-dimensional domains, it has been shown that the compartment size cannot be chosen arbitrarily small for systems containing second-order or higher-order reactions, i.e., the error of a compartment-based simulation increases as the compartment size approaches zero (Erban and Chapman 2020). There have been numerous studies (Isaacson 2009; Erban and Chapman 2009; Hellander et al. 2012; Kang et al. 2012; Isaacson 2013; Agbanusi and Isaacson 2014; Hellander and Petzold 2016; Isaacson and Zhang 2018) on the optimal choice of the compartment size and its influence on the approximation quality compared to microscopic modeling approaches given by the Doi model (Doi 1976a, b) or the Smoluchowski model (Andrews et al. 2010; Andrews 2012).

The optimal size of the compartment depends on the diffusion constant (Erban and Chapman 2009; Isaacson 2009; Hellander et al. 2012; Kang et al. 2012; Isaacson 2013; Agbanusi and Isaacson 2014; Hellander and Petzold 2016; Isaacson and Zhang 2018). In particular, if a biological system consists of molecular species with different diffusion constants, the compartment-based model can be naturally generalized to allow for different meshes (compartment sizes) for different chemical species (Cao and Erban 2014; Hellander and Hellander 2020; Li and Cao 2012). In this paper, we analyze such models using a generalization of the RDME which we will call the *multi-grid reaction diffusion master equation* (mgRDME). Numerical simulations of multi-grid reaction-diffusion models allow for high accuracy at reduced numerical cost compared to fully microscopically resolved systems (Cao and Erban 2014; Chatterjee et al. 2004; Dai et al. 2008). In some computational frameworks (Hellander and Hellander 2020), molecules can transfer from a fine-grained mesh to a coarse-grained mesh whenever appropriate. In the mgRDME that will be analyzed in this work, the grid size remains constant for each species but varies across species.

The mgRDME framework will be analyzed by applying it to biochemical systems of morphogen gradient formation (Saunders and Howard 2009; Kicheva et al. 2007; Bergmann et al. 2007; Gregor et al. 2007; Iyer et al. 2023). Morphogens are signaling molecules whose non-uniform distribution control pattern formation (or morphogenesis) during the development of multicellular organisms. Starting with the pioneering work of Turing (Turing 1952), it has been shown that pattern formation can naturally arise in reaction-diffusion systems through diffusion-driven instability, provided that the system contains at least two chemical species with different diffusion constants (Murray 2002). In particular, stochastic modelling of Turing patterns is a natural application area of the mgRDME framework (Cao and Erban 2014).

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 3 of 31

6

In this paper, we will focus on pattern formation (morphogen gradient formation) which results from pre-patterning, i.e., we assume that the studied domain has already been differentiated into two regions and the release of signalling molecules is localized in one of those regions (Rogers and Schier 2011; Shimmi and O'Connor 2003; Shimmi et al. 2005). The formation of morphogen gradients then leads to further pattern formation as the cells recognize and interpret the high or low morphogen concentration (Wolpert et al. 2015). Stochastic models helps us to understand the impact of noise (fluctuations) on the pattern (Teimouri and Kolomeisky 2018; Kolomeisky 2011; Erban and Chapman 2020).

In our model, we distinguish between signaling molecules \( A \), which are produced locally in a specific spatial region, and morphogen molecules \( B \), which are produced from \( A \) through a reaction that can occur anywhere in space. Both species move in space via diffusion. The locally restricted production of \( A \) (pre-patterning) leads to creation of a morphogen gradient which is then interpreted by cells to produce further patterning (Erban and Chapman 2020, Section 6.8). Since the signaling molecules \( A \) are assumed to diffuse at a significantly higher rate compared to \( B \), our reaction-diffusion systems are well suited to be analyzed via the mgRDME. Moreover, the primary focus lies on the spatial arrangement of the morphogen \( B \), represented by its gradient, while the specific spatial arrangement of signaling molecules \( A \) is of less importance. This prioritization is reflected in the use of distinct compartment sizes for the two species, with a finer resolution dedicated to the morphogen \( B \).

We study two reaction networks. In Sect. 2, we consider a reaction-diffusion system with the two chemical species A (signal) and B (morphogen) being subject to the following three first-order chemical reactions

\[
\emptyset \longrightarrow A \quad (\text { locally   restricted }) \quad \text { and } \quad A \longrightarrow B \longrightarrow \emptyset , \tag {1}
\]

where the empty set  \( \varnothing \)  is interpreted as sources and sinks of molecules (Erban and Chapman 2020), i.e., molecules of signal A are continuously produced and converted into molecules of morphogen B, which are degraded over time. In Sect. 3, we replace the first-order conversion reaction  \( A \rightarrow B \)  by the second-order dimerization reaction  \( A + A \rightarrow B \) , i.e., the chemical reactions (1) are changed to

\[
\emptyset \longrightarrow A \quad (\text { locally   restricted }) \quad \text { and } \quad A + A \longrightarrow B \longrightarrow \emptyset , \tag {2}
\]

while the diffusion part of the reaction-diffusion model remains the same in both Sects. 2 and 3. Since morphogen gradient systems are effectively one-dimensional (with one 'important' direction), we study both systems in a one-dimensional domain. For reaction-diffusion systems with first-order chemical kinetics as in (1), the one-dimensional results are directly applicable in higher dimensions. However, since our second system (2) includes the second-order reaction (dimerization), the diffusion-limited results will depend on the dimension of the physical space (Kang and Redner 1984; Montroll 1969; Ben-Avraham and Redner 1986).

For each of the two reaction networks (1) and (2), we formulate the standard RDME, the mgRDME, as well as a ‘ground truth’ model given by particle-based Brownian dynamics. The compartment-based modeling approaches are compared to

Springer

6

Page 4 of 31

R. Erban, S. Winkelmann

(a)

![img-1.jpeg](img-1.jpeg)

(b)

![img-2.jpeg](img-2.jpeg)

Fig. 1 Standard compartment-based model. a Schematic of the computational domain \([0, L]\) divided into \(K\) compartments of size \(h\). In the first compartment on the left, highlighted in red, molecules of type \(A\) are produced at rate \(k_{1} > 0\). b The average number of molecules of \(A\) (red) and \(B\) (blue) at steady state obtained by solving Eq. (16) for parameters given in (18) and for \(L = 2\) and \(K = 40\) (i.e., we have \(h = L / K = 0.05\), \(D_A / h^2 = 64\) and \(D_B / h^2 = 4\)) (Color figure online)

the Brownian dynamics results based on the differences in their steady-state distributions. Specifically, we define an error function in terms of the distance between the distribution of morphogen molecules B at steady state, and a cost function as a measure of the numerical complexity when simulating the corresponding stochastic processes. Minimizing error and cost over possible choices of compartment sizes for A and B gives rise to a multi-objective optimization problem, which will be studied for both reaction networks. Our results are further summarized in the discussion Sect. 4. The time-dependent problems are discussed in Appendices 5.1 and 5.4.

## 2 First-Order Fast-Slow Morphogen Gradient Model

In this section, we consider a reaction-diffusion system with the first-order chemical reactions (1). We derive analytic expressions for the long-term spatial distribution of particles, examining the standard RDME, the mgRDME and the Brownian dynamics models in Sects. 2.1, 2.2 and 2.3, respectively. Finally, we compare the models in Sect. 2.4 by means of a multi-objective optimization problem.

### 2.1 Standard Compartment-Based Model

We consider molecules of chemical species A and B that diffuse in the domain  \( [0, L] \) , where L is the domain length. The diffusion coefficients are denoted  \( D_{A} \)  and  \( D_{B} \) ,

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 5 of 31

6

respectively. They have physical units of  \( [length]^{2}/[time] \) . As for the standard RDME discretization, we divide the domain  \( [0, L] \)  into K = L/h compartments of equal length h > 0. This domain is schematically shown in Fig. 1a. Denoting the chemical species A (resp. B) in the i-th compartment  \( [(i - 1)h, ih) \)  by  \( A_{i} \)  (resp.  \( B_{i} \) ), where  \( i = 1, 2, \ldots, K \) , then diffusion corresponds to two chains of "chemical reactions" (Cao and Erban 2014; Erban and Chapman 2020):

\[
A _ {1} \xrightarrow [ d _ {A} ]{d _ {A}} A _ {2} \xrightarrow [ d _ {A} ]{d _ {A}} A _ {3} \xrightarrow [ d _ {A} ]{d _ {A}} \dots \xrightarrow [ d _ {A} ]{d _ {A}} A _ {K}, \tag {3}
\]

\[
B _ {1} \xrightarrow [ d _ {B} ]{d _ {B}} B _ {2} \xrightarrow [ d _ {B} ]{d _ {B}} B _ {3} \xrightarrow [ d _ {B} ]{d _ {B}} \dots \xrightarrow [ d _ {B} ]{d _ {B}} B _ {K}, \tag {4}
\]

where

\[
d _ {A} = \frac {D _ {A}}{h ^ {2}} \quad \text { and } \quad d _ {B} = \frac {D _ {B}}{h ^ {2}} \tag {5}
\]

are the jump rates. Since h has the physical unit of [length], the jump rates have units  \( [time]^{-1} \) . We assume that molecules of species A are produced in the first compartment on the left, highlighted in red in Fig. 1a. This corresponds to the reaction

\[
\emptyset \xrightarrow {k _ {1} / h} A _ {1}, \tag {6}
\]

where the rate constant  \( k_{1} > 0 \)  has units  \( [time]^{-1} \) . Scaling with h in reaction (6) is necessary to ensure that after multiplying by the compartment size h – as is usual for zero-order reactions (Erban and Chapman 2020) – we get a production rate independent of the compartment size. This gives the same production rate in the first compartment for any value of h, and is consistent with our ‘ground truth’ Brownian dynamics model corresponding to partial differential equations (PDEs) (27)–(28) in Sect. 2.3. Since (1) assumes that A is converted to B, and B is degraded in the whole domain, we have the reactions

\[
A _ {i} \xrightarrow {k _ {2}} B _ {i} \xrightarrow {k _ {3}} \emptyset , \quad \text { for } i = 1, 2, \dots , K, \tag {7}
\]

with rate constants  \( k_{2} \)  and  \( k_{3} \)  again having physical units of  \( [time]^{-1} \) . The (random) number of particles in the i-th compartment at time  \( t \geq 0 \)  is denoted by  \( A_{i}(t) \)  and  \( B_{i}(t) \) , respectively. Let  \( p(\mathbf{n}, \mathbf{m}, t) \)  be the joint probability that  \( A_{i}(t) = n_{i} \)  and  \( B_{i}(t) = m_{i} \)  for  \( i = 1, 2, \ldots, K \) , where we use the notation  \( n = [n_{1}, n_{2}, \ldots, n_{K}] \in N^{K} \)  and  \( m = [m_{1}, m_{2}, \ldots, m_{K}] \in N^{K} \)  to define the system state. To formulate the RDME corresponding to the chemical reaction system (3)–(7), we define operators  \( O_{i}^{+}, O_{i}^{-}: N^{K} \to N^{K} \)  by

\[
\mathcal {O} _ {i} ^ {+} \left[ n _ {1}, \dots , n _ {i - 1}, n _ {i}, n _ {i + 1}, \dots , n _ {K} \right] := \left[ n _ {1}, \dots , n _ {i - 1}, n _ {i} + 1, n _ {i + 1}, \dots , n _ {K} \right], \tag {8}
\]

\[
\mathcal {O} _ {i} ^ {-} \left[ n _ {1}, \dots , n _ {i - 1}, n _ {i}, n _ {i + 1}, \dots , n _ {K} \right] := \left[ n _ {1}, \dots , n _ {i - 1}, n _ {i} - 1, n _ {i + 1}, \dots , n _ {K} \right], \tag {9}
\]

Springer

6

Page 6 of 31

R. Erban, S. Winkelmann

for \(i = 1,2,\ldots ,K\). By means of these operators, we define the diffusion operator \(\mathcal{D}:L^{1}\big(\mathbb{N}^{K}\times \mathbb{N}^{K}\big)\to L^{1}\big(\mathbb{N}^{K}\times \mathbb{N}^{K}\big)\), where

\[
L ^ {1} \left(\mathbb {N} ^ {K} \times \mathbb {N} ^ {K}\right) := \left\{f: \mathbb {N} ^ {K} \times \mathbb {N} ^ {K} \rightarrow \mathbb {R} \mid \sum_ {\mathbf {n}, \mathbf {m}} f (\mathbf {n}, \mathbf {m}) <   \infty \right\} \tag {10}
\]

by

\[
\begin{array}{l} \mathcal {D} f (\mathbf {n}, \mathbf {m}) := \frac {D _ {A}}{h ^ {2}} \sum_ {i = 1} ^ {K - 1} \left\{(n _ {i} + 1) f (\mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i + 1} ^ {-} \mathbf {n}, \mathbf {m}) - n _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {A}}{h ^ {2}} \sum_ {i = 2} ^ {K} \left\{(n _ {i} + 1) f (\mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i - 1} ^ {-} \mathbf {n}, \mathbf {m}) - n _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {B}}{h ^ {2}} \sum_ {i = 1} ^ {K - 1} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i + 1} ^ {-} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {B}}{h ^ {2}} \sum_ {i = 2} ^ {K} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i - 1} ^ {-} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \tag {11} \\ \end{array}
\]

and the reaction operator \(\mathcal{R}:L^{1}\big(\mathbb{N}^{K}\times \mathbb{N}^{K}\big)\to L^{1}\big(\mathbb{N}^{K}\times \mathbb{N}^{K}\big)\) by

\[
\begin{array}{l} \mathcal {R} f (\mathbf {n}, \mathbf {m}) := k _ {1} \left\{f \left(\mathcal {O} _ {1} ^ {-} \mathbf {n}, \mathbf {m}\right) - f (\mathbf {n}, \mathbf {m}) \right\} \\ + k _ {2} \sum_ {i = 1} ^ {K} \left\{(n _ {i} + 1) f (\mathcal {O} _ {i} ^ {+} \mathbf {n}, \mathcal {O} _ {i} ^ {-} \mathbf {m}) - n _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ + k _ {3} \sum_ {i = 1} ^ {K} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\}. \tag {12} \\ \end{array}
\]

Then the reaction-diffusion master equation, which corresponds to the system of reactions (3-7), can be written as follows

\[
\frac {\partial p}{\partial t} (\mathbf {n}, \mathbf {m}, t) = (\mathcal {D} + \mathcal {R})   p (\mathbf {n}, \mathbf {m}, t). \tag {13}
\]

The stationary distribution is defined by

\[
\phi (\mathbf {n}, \mathbf {m}) = \lim _ {t \rightarrow \infty} p (\mathbf {n}, \mathbf {m}, t)
\]

and satisfies the stationary reaction-diffusion master equation

\[
0 = (\mathcal {D} + \mathcal {R}) \phi (\mathbf {n}, \mathbf {m}). \tag {14}
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 7 of 31

6

Solving (14), we obtain the product of Poisson distributions (Jahnke and Huisinga 2007)

\[
\phi (\mathbf {n}, \mathbf {m}) = \exp \left[ - \sum_ {i = 1} ^ {K} \left(\overline {{A}} _ {i} + \overline {{B}} _ {i}\right) \right] \prod_ {i = 1} ^ {K} \frac {\overline {{A}} _ {i} ^ {n _ {i}} \overline {{B}} _ {i} ^ {m _ {i}}}{n _ {i} ! m _ {i} !}, \tag {15}
\]

where  \( \overline{A}_{i} \)  and  \( \overline{B}_{i} \) , for  \( i = 1, 2, \ldots, K \) , satisfy the steady-state equations

\[
\left(\frac {D _ {A}}{h ^ {2}} S - k _ {2} I\right) \overline {{{\mathbf {A}}}} = - k _ {1} \mathbf {e} _ {1}, \quad \left(\frac {D _ {B}}{h ^ {2}} S - k _ {3} I\right) \overline {{{\mathbf {B}}}} = - k _ {2} \overline {{{\mathbf {A}}}}. \tag {16}
\]

Here, \(I \in \mathbb{R}^{K \times K}\) is the identity matrix and vectors \(\overline{\mathbf{A}} \in \mathbb{R}^K\), \(\overline{\mathbf{B}} \in \mathbb{R}^K\), \(\mathbf{e}_1 \in \mathbb{R}^K\) and matrix \(S \in \mathbb{R}^{K \times K}\) are given by

\[
\overline {{{\mathbf {A}}}} = \left( \begin{array}{c} \overline {{{A}}} _ {1} \\ \overline {{{A}}} _ {2} \\ \overline {{{A}}} _ {3} \\ \vdots \\ \overline {{{A}}} _ {K - 1} \\ \overline {{{A}}} _ {K} \end{array} \right), \quad \overline {{{\mathbf {B}}}} = \left( \begin{array}{c} \overline {{{B}}} _ {1} \\ \overline {{{B}}} _ {2} \\ \overline {{{B}}} _ {3} \\ \vdots \\ \overline {{{B}}} _ {K - 1} \\ \overline {{{B}}} _ {K} \end{array} \right), \quad \mathbf {e} _ {1} = \left( \begin{array}{c} 1 \\ 0 \\ 0 \\ \vdots \\ 0 \\ 0 \end{array} \right), \quad S = \left( \begin{array}{c c c c c c c} - 1 & 1 & 0 & 0 & \dots & 0 & 0 \\ 1 & - 2 & 1 & 0 & \dots & 0 & 0 \\ 0 & 1 & - 2 & 1 & \dots & 0 & 0 \\ \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\ 0 & 0 & 0 & 0 & \dots & - 2 & 1 \\ 0 & 0 & 0 & 0 & \dots & 1 & - 1 \end{array} \right). \tag {17}
\]

The entries  \( \overline{A}_{i} \)  and  \( \overline{B}_{i} \)  are the average numbers of molecules of A and B, respectively, in the i-th compartment at equilibrium. The solution  \( (\overline{\mathbf{A}},\overline{\mathbf{B}}) \)  of Eq. (16) is plotted in Fig. 1b for the following choice of rate parameter values

\[
k _ {1} = 1 0 0, \quad k _ {2} = 2, \quad k _ {3} = 1, \quad D _ {A} = 0. 1 6, \quad D _ {B} = 0. 0 1, \tag {18}
\]

showing the morphogen gradient. As the stationary distribution  \( \phi \)  is given by a product of Poisson distributions in (15), the values  \( \overline{A}_{i} \)  and  \( \overline{B}_{i} \)  not only determine the long-term averages but also the variances of the particle numbers at steady state. Since the values of  \( \overline{A}_{i} \)  and  \( \overline{B}_{i} \)  fully characterize the stationary distribution  \( \phi \) , we will compare the models in terms of  \( \overline{A}_{i} \)  and  \( \overline{B}_{i} \) .

### 2.2 Generalized Compartment-Based Model: mgRDME

Assuming that molecules of species A diffuse significantly faster than those of species B (as indicated by our parameter values in Eq. (18), where  \( D_{A} \gg D_{B} \) ), we choose a larger compartment size for species A (Cao and Erban 2014), while keeping the compartments for species B of the same size as in the standard model in Sect. 2.1. Let  \( h_{A} = L/K_{A} \)  for  \( K_{A} \in N \)  be the compartment length for species A, and  \( h_{B} = L/K_{B} \)  for  \( K_{B} \in N \) , where  \( K_{B} > K_{A} \) , be the compartment length for species B. For simplicity, we will assume that the ratio  \( \gamma := K_{B}/K_{A} = h_{A}/h_{B} > 1 \)  is a natural number, i.e.  \( \gamma \in N \) . However, a generalization is possible, see Remark 1 for an explanation. In the mgRDME formulation (Cao and Erban 2014), reaction (6) is replaced by

Springer

6

Page 8 of 31

R. Erban, S. Winkelmann

\[
\emptyset \xrightarrow {k _ {1} / h _ {A}} A _ {1}, \tag {19}
\]

which implies again that \( A \)-particles are produced at rate \( k_{1} \) in the first compartment on the left (of the coarser \( A \)-discretization), independently of the grid size \( h_{A} \). The reaction \( A_{i} \longrightarrow B_{i} \) in (7) is generalized to

\[
A _ {j} \xrightarrow {k _ {2} / \gamma} B _ {i}, \quad j = 1, 2, \dots , K _ {A}, \quad \text { for } i \in \mathcal {I} (j) := \{(j - 1) \gamma + 1, \dots , j \gamma \}, \tag {20}
\]

which means that after conversion from A to B the resulting B-molecule is placed uniformly in one of the  \( \gamma \)  (smaller) compartments that overlap with the j-th (large) compartment of the reacting A-particle. The diffusion operator (11) is generalized to

\[
\begin{array}{l} \tilde {\mathcal {D}} f (\mathbf {n}, \mathbf {m}) := \frac {D _ {A}}{h _ {A} ^ {2}} \sum_ {j = 1} ^ {K _ {A} - 1} \left\{(n _ {j} + 1) f (\mathcal {O} _ {j} ^ {+} \mathcal {O} _ {j + 1} ^ {-} \mathbf {n}, \mathbf {m}) - n _ {j} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {A}}{h _ {A} ^ {2}} \sum_ {j = 2} ^ {K _ {A}} \left\{(n _ {j} + 1) f (\mathcal {O} _ {j} ^ {+} \mathcal {O} _ {j - 1} ^ {-} \mathbf {n}, \mathbf {m}) - n _ {j} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {B}}{h _ {B} ^ {2}} \sum_ {i = 1} ^ {K _ {B} - 1} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i + 1} ^ {-} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {D _ {B}}{h _ {B} ^ {2}} \sum_ {i = 2} ^ {K _ {B}} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i - 1} ^ {-} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\}, \tag {21} \\ \end{array}
\]

while the reaction operator (12) now reads

\[
\begin{array}{l} \tilde {\mathcal {R}} f (\mathbf {n}, \mathbf {m}) := k _ {1} \left\{f \left(\mathcal {O} _ {1} ^ {-} \mathbf {n}, \mathbf {m}\right) - f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {k _ {2}}{\gamma} \sum_ {j = 1} ^ {K _ {A}} \sum_ {i \in \mathcal {I} (j)} \left\{(n _ {j} + 1) f (\mathcal {O} _ {j} ^ {+} \mathbf {n}, \mathcal {O} _ {i} ^ {-} \mathbf {m}) - n _ {j} f (\mathbf {n}, \mathbf {m}) \right\} \\ + k _ {3} \sum_ {i = 1} ^ {K _ {B}} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\}. \tag {22} \\ \end{array}
\]

Dividing by  \( \gamma \)  in the second line of (22) goes in line with summing over the set  \( \mathcal{I}(j) \)  defined in (20) which contains  \( |\mathcal{I}(j)| = \gamma \)  elements; this can be interpreted as a sum of  \( \gamma \)  reactions  \( A_{j} \to B_{i} \)  that share the rate  \( k_{2} \) . In total, we obtain the multi-grid reaction-diffusion master equation (mgRDME)

\[
\frac {\partial p}{\partial t} (\mathbf {n}, \mathbf {m}, t) = (\tilde {\mathcal {D}} + \tilde {\mathcal {R}}) p (\mathbf {n}, \mathbf {m}, t). \tag {23}
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 9 of 31

6

Solving the corresponding stationary mgRDME, given by  \( 0 = (\tilde{\mathcal{D}} + \tilde{\mathcal{R}}) \phi(\mathbf{n}, \mathbf{m}) \)  in analogy to (14), we obtain a product of Poisson distributions similar to (15):

\[
\phi (\mathbf {n}, \mathbf {m}) = \exp \left[ - \sum_ {j = 1} ^ {K _ {A}} \overline {{A}} _ {j} - \sum_ {i = 1} ^ {K _ {B}} \overline {{B}} _ {i} \right] \prod_ {j = 1} ^ {K _ {A}} \frac {\overline {{A}} _ {j} ^ {n _ {j}}}{n _ {j} !} \prod_ {i = 1} ^ {K _ {B}} \frac {\overline {{B}} _ {i} ^ {m _ {i}}}{m _ {i} !},
\]

where \(\overline{A}_j\) and \(\overline{B}_i\), for \(j = 1,2,\ldots,K_A\), \(i = 1,2,\ldots,K_B\), satisfy the generalized steady-state equations

\[
\left(\frac {D _ {A}}{h _ {A} ^ {2}} S _ {A} - k _ {2} I _ {A}\right) \overline {{{\mathbf {A}}}} = - k _ {1} \mathbf {e} _ {1}, \quad \left(\frac {D _ {B}}{h _ {B} ^ {2}} S _ {B} - k _ {3} I _ {B}\right) \overline {{{\mathbf {B}}}} = - \frac {k _ {2}}{\gamma} M \overline {{{\mathbf {A}}}}. \tag {24}
\]

Here, \( I_A \in \mathbb{R}^{K_A \times K_A} \) and \( I_B \in \mathbb{R}^{K_B \times K_B} \) are identity matrices and vectors \( \overline{\mathbf{A}} \in \mathbb{R}^{K_A} \), \( \overline{\mathbf{B}} \in \mathbb{R}^{K_B} \), \( \mathbf{e}_1 \in \mathbb{R}^{K_A} \) and matrices \( S_A \in \mathbb{R}^{K_A \times K_A} \), \( S_B \in \mathbb{R}^{K_B \times K_B} \) are defined in an analogous manner as before in (17). In addition, there is the block matrix \( M \in \mathbb{R}^{K_B \times K_A} \) given by

\[
M = \left( \begin{array}{c} M _ {1} \\ \vdots \\ M _ {K _ {A}} \end{array} \right), \quad M _ {j} = \left( \begin{array}{c c c c c c} 0 & \ldots & 0 & 1 & 0 & \ldots & 0 \\ \vdots & & \vdots & \vdots & \vdots & & \vdots \\ 0 & \ldots & 0 & 1 & 0 & \ldots & 0 \end{array} \right) \in \mathbb {R} ^ {\gamma \times K _ {A}}, \quad j = 1, \ldots , K _ {A},
\]

with the non-zero entries (i.e., the ones) in the block  \( M_{j} \in R^{\gamma \times K_{A}} \)  placed in the j-th column. As in the standard setting,  \( \overline{A}_{j} \)  and  \( \overline{B}_{i} \)  give both the long-term mean and variance of the population sizes in the respective boxes, and they fully characterize the stationary distribution  \( \phi \) .

Remark 1 The choice of an integer-valued ratio \(\gamma\) simplifies notation in (20) and (22). However, a generalization to an arbitrary ratio \(K_B / K_A \notin \mathbb{N}\) is straightforward: Let \(C_j^{(A)} := [(j-1)h_A, jh_A)\) and \(C_i^{(B)} := [(i-1)h_B, ih_B)\) denote the \(j\)-th and \(i\)-th compartment of the discretizations chosen for species \(A\) and species \(B\), respectively. Then, reaction \(A_j \to B_i\), meaning that a particle of species \(A\) located in \(C_j^{(A)}\) is converted into \(B\) located in \(C_i^{(B)}\), occurs at rate

\[
\frac {k _ {2} \left| C _ {i} ^ {(B)} \cap C _ {j} ^ {(A)} \right|}{\left| C _ {j} ^ {(A)} \right|}, \tag {25}
\]

where  \( |\cdot| \)  represents to the size of the compartment, specifically the length of the interval in the one-dimensional problems considered here. This means that the conversion rate  \( k_{2} \)  is shared among the potential target compartments according to their overlap with the source compartment.

Springer

6

Page 10 of 31

R. Erban, S. Winkelmann

![img-3.jpeg](img-3.jpeg)

![img-4.jpeg](img-4.jpeg)

(a) \(\gamma = 2\)

![img-5.jpeg](img-5.jpeg)

![img-6.jpeg](img-6.jpeg)

(b) \(\gamma = 4\)

Fig. 2 Standard RDME versus mgRDME. Average number of molecules of A and B at equilibrium obtained by solving the standard stationary RDME (16) (red/blue dots) and by solving the stationary mgRDME (24) (grey bars) for rate constants given in (18) and for L = 2,  \( K_{B} = 40 \)  and a  \( \gamma = 2 \) , b  \( \gamma = 4 \)  (Color figure online)

![img-7.jpeg](img-7.jpeg)

Fig. 3 Difference between standard RDME and mgRDME. Difference \(\| \bar{\mathbf{B}}_1 - \bar{\mathbf{B}}_\gamma \|_1\) as a function of the ratio \(\gamma = K_B / K_A\) for rate constants given in (18) and for \(L = 2\), \(K_B = 60\) and \(K_A \in \{60, 30, 20, 15, 12, 10, 6, 4, 3\}\), where \(\bar{\mathbf{B}}_1\) the solution of the standard steady-state Eq. (16), and \(\bar{\mathbf{B}}_\gamma\) the solution of the generalized steady-state Eq. (24) (Color figure online)

#### Model comparison: standard RDME versus mgRDME.

In the following, we denote the solution of the generalized steady-state Eq. (24) by \(\bar{\mathbf{A}}_{\gamma} = (\bar{A}_j^{(\gamma)})_{j=1,2,\dots,K_A}\), \(\bar{\mathbf{B}}_{\gamma} = (\bar{B}_i^{(\gamma)})_{i=1,2,\dots,K_B}\) in order to emphasize its dependence on the ratio \(\gamma = K_B / K_A\) and to distinguish from the solution \(\bar{\mathbf{A}} = \bar{\mathbf{A}}_1\), \(\bar{\mathbf{B}} = \bar{\mathbf{B}}_1\) of the standard compartment-based model given by Eq. (16). Figure 2 shows the solution \(\bar{\mathbf{A}}_{\gamma}\), \(\bar{\mathbf{B}}_{\gamma}\) of equation (24) for the same parameter values as in Fig. 1b and for \(\gamma = 2\) and \(\gamma = 4\). Instead of \(\bar{\mathbf{A}}_{\gamma}\) we plot the rescaled values \(\bar{A}_i^{(\gamma)} := \bar{A}_j^{(\gamma)} / \gamma\) for \(i \in \mathcal{I}(j)\). This is for the purpose of comparability with the solution \(\bar{\mathbf{A}}\) of the standard steady-state Eq. (16) plotted in Fig. 1b.

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 11 of 31

6

We observe a qualitatively good agreement of the generalized solutions with the original one. To see how this agreement depends on the level of coarsening, we plot the difference  \( \|\bar{B}_{1}-\bar{B}_{\gamma}\|_{1} \)  as a function of  \( \gamma=K_{B}/K_{A} \)  in Fig. 3. We choose  \( \|\bar{B}_{1}-\bar{B}_{\gamma}\|_{1} \)  as the measure of error because the population of species B (morphogen) is the output of interest. Note that  \( \bar{B}_{\gamma} \)  is a vector of fixed size  \( K_{B} \)  for any  \( \gamma \), whereas the dimension of  \( \bar{A}_{\gamma} \)  is  \( K_{A} \), which depends on  \( \gamma \). Consequently, the distance between  \( \bar{A}_{1} \)  and  \( \bar{A}_{\gamma} \)  is not directly defined for  \( \gamma\neq1 \). We see from Fig. 3 that the difference  \( \|\bar{B}_{1}-\bar{B}_{\gamma}\|_{1} \)  increases monotonically with  \( \gamma \), which is to be expected since larger values of  \( \gamma \)  represent a higher degree of coarse-graining. However, the increase of the error is nonlinear, and there is an area of values ( \( \gamma\lesssim5 \) ) where the error is sufficiently small to be deemed acceptable. Before we study the approximation quality in more detail, we define in the following section the ‘ground truth’ model given by particle-based Brownian dynamics.

### 2.3 'Ground Truth' Particle-Based Model: Brownian Dynamics

As our ‘ground truth’ model, we will consider a model with two species A and B diffusing and reacting in the unbounded domain R (Erban and Chapman 2020). Molecules of A are released at the origin, x = 0, with rate  \( 2k_{1} \) , which has units  \( [sec^{-1}] \) . Due to symmetry, this is equivalent to studying the same model on the half-line  \( [0, \infty) \)  with the rate constant  \( k_{1} \)  and reflecting boundary conditions at x = 0 (Erban and Chapman 2007). Molecules of A and B are diffusing with diffusion constants  \( D_{A} \)  and  \( D_{B} \) , respectively. They are subject to the conversion and degradation reactions

\[
A \xrightarrow {k _ {2}} B \xrightarrow {k _ {3}} \emptyset . \tag {26}
\]

We compare the Brownian dynamics to the simulation on the finite interval  \( [0, L] \) , which is discretized into compartments in the RDME and mgRDME frameworks. The effect of the boundary at x = L is discussed in Sect. 5.2 of the Appendix.

Let \( a(x, t) \) and \( b(x, t) \) be the concentrations of molecules of \( A \) and \( B \), respectively, that is \( a(x, t) \, \mathrm{d}x \) denotes the average number of molecules in the interval \( [x, x + \mathrm{d}x) \) at time \( t \) in the Brownian dynamics model. Then \( a(x, t) \) and \( b(x, t) \) satisfy the following PDEs (Kicheva et al. 2007, equation (1))

\[
\frac {\partial a}{\partial t} = D _ {A} \frac {\partial^ {2} a}{\partial x ^ {2}} - k _ {2} a + 2 k _ {1} \delta (x), \tag {27}
\]

\[
\frac {\partial b}{\partial t} = D _ {B} \frac {\partial^ {2} b}{\partial x ^ {2}} + k _ {2} a - k _ {3} b, \tag {28}
\]

where \(x \in \mathbb{R}\), and \(\delta(x)\) is the Dirac delta function. We consider the boundary conditions

\[
\lim _ {x \to \pm \infty} a (x, t) = \lim _ {x \to \pm \infty} b (x, t) = 0. \tag {29}
\]

Springer

6

Page 12 of 31

R. Erban, S. Winkelmann

Let \(\bar{a}(x)\) and \(\bar{b}(x)\) be the corresponding stationary distributions, i.e.

\[
\bar {a} (x) = \lim _ {t \to \infty} a (x, t), \quad \bar {b} (x) = \lim _ {t \to \infty} b (x, t).
\]

They satisfy the steady-state equations

\[
0 = D _ {A} \frac {\mathrm{d} ^ {2} \bar {a}}{\mathrm{d} x ^ {2}} - k _ {2} \bar {a} + 2 k _ {1} \delta (x), \tag {30}
\]

\[
0 = D _ {B} \frac {\mathrm{d} ^ {2} \bar {b}}{\mathrm{d} x ^ {2}} + k _ {2} \bar {a} - k _ {3} \bar {b}, \tag {31}
\]

with boundary conditions

\[
\lim _ {x \to \pm \infty} \bar {a} (x) = \lim _ {x \to \pm \infty} \bar {b} (x) = 0. \tag {32}
\]

We notice that the PDE (27) (resp. the ODE (30)) does not depend on  \( b(x, t) \)  (resp.  \( \bar{b}(x) \) ). In what follows, we will therefore first solve the Eq. (27) for the chemical species A, and then we will substitute it into the second Eq. (28) for the chemical species B to get our ‘ground truth’ solution. In the following, we will do this for the steady-state problem. The time-dependent problem is considered in the Appendix, see Sect. 5.1.

#### Steady-state solution

Since the second derivative of the absolute value function, \( |x| \), is equal to \( 2\delta(x) \), we can write the solution to Eq. (30) with boundary conditions (32) as

\[
\bar {a} (x) = \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | x | \right]. \tag {33}
\]

This result is also obtained when taking the limit \( t \to \infty \) of the time-dependent solution \( a(x, t) \) as given in Eq. (49) of Appendix 5.1. Substituting into Eq. (31), we get

\[
D _ {B} \frac {\mathrm{d} ^ {2} \bar {b}}{\mathrm{d} x ^ {2}} + k _ {1} \sqrt {\frac {k _ {2}}{D _ {A}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | x | \right] - k _ {3} \bar {b} = 0.
\]

Consequently, using boundary conditions (32), we obtain

\[
\begin{array}{l} \bar {b} (x) = \frac {k _ {1} \sqrt {k _ {2} D _ {A}}}{D _ {A} k _ {3} - D _ {B} k _ {2}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | x | \right] \\ + \frac {k _ {1} k _ {2} \sqrt {D _ {B}}}{\left(D _ {B} k _ {2} - D _ {A} k _ {3}\right) \sqrt {k _ {3}}} \exp \left[ - \sqrt {\frac {k _ {3}}{D _ {B}}} | x | \right] \tag {34} \\ \end{array}
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 13 of 31

6

for the concentration  \( \bar{b} \)  of molecules of B at equilibrium. Our compartment-based models simulate the problem on the finite interval  \( [0, L] \) , with zero-flux boundary condition at x = L. To get the ‘ground truth’ solution (34), we have used the infinite domain, R, with boundary conditions (32). In Sect. 5.2 we show that the effect of non-flux boundary conditions on the ‘ground truth’ solution  \( \bar{a} \)  vanishes for  \( L \to \infty \) . In particular, the concentration is negligible at  \( x \geq L \)  as long as L is large enough, and we can use (33) and (34) for a comparison with the compartment-based models.

### 2.4 Model Comparison

Since the reaction system under consideration is a first-order reaction network, the standard compartment-based model of Sect. 2.1 converges to the Brownian dynamics solution for  \( h \rightarrow 0 \)  (Engblom et al. 2009; Erban and Chapman 2020). Our goal is now to compare Brownian dynamics with the compartment-based model for finite h, both for the standard RDME and for the mgRDME. In application, one is mainly interested in the output given by the product B, which motivates to define the distance between the models by means of the average number of molecules of B.

#### 2.4.1 Error Definition

Given the steady-state solution  \( \bar{\mathbf{B}} = (\bar{B}_{1}, \bar{B}_{2}, \ldots, \bar{B}_{K_{B}}) \)  of the mgRDME model with  \( K_{B} \)  boxes of size  \( h = L/K_{B} \)  (which agrees with the standard RDME solution when choosing  \( \gamma = 1 \) ), we define the piecewise constant function  \( \bar{b}_{h} : [0, L] \to [0, \infty) \) ,

\[
\bar {b} _ {h} (x) := \bar {B} _ {i} / h \quad \text { for } x \in [ (i - 1) h, i h), \quad i = 1, 2, \dots , K _ {B}, \tag {35}
\]

which approximates the steady-state concentration of molecules of B as it depends on the location x. The  \( L^{1} \) -distance to the Brownian dynamics solution  \( \bar{b} \)  defines

\[
\operatorname{err} \left(K _ {B}, \gamma\right) := \left| \left| \bar {b} - \bar {b} _ {h} \right| \right| _ {L ^ {1}} = \int_ {0} ^ {L} | \bar {b} (x) - \bar {b} _ {h} (x) | \mathrm{d} x = \sum_ {i = 1} ^ {K _ {B}} \int_ {(i - 1) h} ^ {i h} | \bar {b} (x) - \bar {B} _ {i} | \mathrm{d} x, \tag {36}
\]

the spatial error in dependence on  \( K_{B} \)  and on the ratio  \( \gamma = K_{B}/K_{A} \) . We note that this error would be positive even if  \( \bar{b}_{h} \)  was a discretization of  \( \bar{b} \)  (and not the solution of the compartment-based model), simply because of the coarse-graining. We thus expect the error to naturally decrease with decreasing h.

Figure 4 shows this error in dependence on the grid size h, where  \( h_{A} = h_{B} = h \)  for the standard RDME, while  \( h_{A} = h \)  and fixed  \( h_{B} = h^{*} = 1/60 \)  for the generalized mgRDME. For both scenarios, there is a monotone increase of the error with h, which is to be expected in a first-order reaction system. The generalized model shows better results as compared to the standard one, which is due to the fixed small grid size  \( h^{*} \)  for B. However, the error  \( \text{err}(K_{B}, \gamma) \)  seems to be linear in h for  \( \gamma = 1 \)  (standard model, blue crosses), while it is nonlinear in  \( h = h_{A} \)  for the generalized system and stays close to zero for small  \( h = h_{A} \) . We can say that for  \( h = h_{A} \leq 0.1 (\gamma \leq 6) \)  the generalized

Springer

6

Page 14 of 31

R. Erban, S. Winkelmann

![img-8.jpeg](img-8.jpeg)

Fig. 4 Error in dependence on the compartment size. Error \(\mathrm{err}(K_B, \gamma)\), defined in (36), between the steady-state solution \(\bar{b}\) of the 'ground truth' model given by (34) and its approximation by the standard and the generalized compartment based models with solutions \(\bar{b}_h\) defined in (35) in dependence on the compartment size \(h\). We have \(h_A = h_B = h\) for the standard model, while \(h_A = h\) and \(h_B = h^* = 1/60\) for the multi-grid model. Rate constants given in (18) and \(L = 2\) (Color figure online)

mgRDME model gives a good approximation, while the standard compartment-based model already shows a clear error.

Despite the availability of an analytical solution for the steady-state distribution, there remains a fundamental interest in stochastic simulations of the dynamics, which allow for a more comprehensive exploration of the system's behavior including complex interactions and rare events. We point out that, in Fig. 4 and for a given \( h \), the numerical effort for simulating trajectories is larger for the multi-grid model compared to the standard model due to the retention of the small grid size \( h_B = h^* \) for \( B \). However, the overall goal of the generalized method is to reduce numerical cost while simultaneously minimizing errors. Therefore, determining an optimal value for \( \gamma \) entails multi-objective optimization, a task we will address next.

#### 2.4.2 Multi-Objective Optimization

In this section, we study the dependence of the spatial error  \( \text{err}(K_{B}, \gamma) \)  defined in (36) and of the numerical effort on the grid sizes  \( h_{A} \)  and  \( h_{B} \) . We calculate error and cost for different parameter values to show the advantage of the generalized method over the standard method. To define the cost function, we use propensities at steady state because these indicate the runtime of our numerical simulations.

##### Propensities at steady state

The total number of particles at steady state is given by, \( ^{1} \)

\[
\bar {A} _ {\text { total }} = \sum_ {i = 1} ^ {K _ {A}} \bar {A} _ {i} = \frac {k _ {1}}{k _ {2}}, \quad \bar {B} _ {\text { total }} = \sum_ {i = 1} ^ {K _ {B}} \bar {B} _ {i} = \frac {k _ {2}}{k _ {3}} \bar {A} _ {\text { total }} = \frac {k _ {1}}{k _ {3}}. \tag {37}
\]

The overall propensity for a reaction to take place at steady state (no matter which reaction or where in space) is

\[
\bar {r} _ {\text { reac }} := k _ {1} + k _ {2} \bar {A} _ {\text { total }} + k _ {3} \bar {B} _ {\text { total }} = 3 k _ {1}. \tag {38}
\]

\( ^{1} \)  The values result from setting the related reaction rate equations to zero, i.e. from solving  \( k_{1}-k_{2}\bar{A}_{total}=0 \)  and  \( k_{2}\bar{A}_{total}-k_{3}\bar{B}_{total}=0 \) . We note that this approach is only applicable to first-order reaction systems.

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 15 of 31

6

Analogously, the overall propensity for a diffusive jump to take place at steady state is

\[
\begin{array}{l} \bar {r} _ {\text { diff }} := \frac {D _ {A}}{h _ {A} ^ {2}} (2 \bar {A} _ {\text { total }} - \bar {A} _ {1} - \bar {A} _ {K _ {A}}) + \frac {D _ {B}}{h _ {B} ^ {2}} (2 \bar {B} _ {\text { total }} - \bar {B} _ {1} - \bar {B} _ {K _ {B}}) \\ \lesssim 2 \frac {D _ {A}}{h _ {A} ^ {2}} \bar {A} _ {\text { total }} + 2 \frac {D _ {B}}{h _ {B} ^ {2}} \bar {B} _ {\text { total }} = 2 \frac {D _ {A}}{h _ {A} ^ {2}} \frac {k _ {1}}{k _ {2}} + 2 \frac {D _ {B}}{h _ {B} ^ {2}} \frac {k _ {1}}{k _ {3}} \tag {39} \\ \end{array}
\]

with the negative terms in the first line resulting from the fact that diffusive jumps at the outer boundaries of the domain can only go in one direction. The latter expression is a good approximation as long as the population sizes in the outer boxes are not too large (as compared to the total population size), which is the case in the scenarios considered.

The number of iterations steps needed to create a trajectory of the dynamics (starting in steady state, using the Gillespie algorithm (Gillespie 1977, 2007)) scales with \(\bar{r}_{\mathrm{reac}} + \bar{r}_{\mathrm{diff}}\). This motivates to define the cost function as a sum of (38) and (39), giving

\[
c (K _ {B}, \gamma) := 3 k _ {1} + 2 \frac {D _ {A}}{h _ {A} ^ {2}} \frac {k _ {1}}{k _ {2}} + 2 \frac {D _ {B}}{h _ {B} ^ {2}} \frac {k _ {1}}{k _ {3}}, \tag {40}
\]

where \( h_B = L / K_B \) and \( h_A = \gamma h_B \). Figure 5 shows the values of error and cost depending on \( K_B \) and \( \gamma \) (where \( \gamma = 1 \) corresponds to the standard RDME model). We see that the generalized model clearly outperforms the standard model: The Pareto front exclusively comprises points derived from the generalized model. Consequently, the decrease in spatial resolution for species \( A \) introduced by the generalized model minimally impacts the error in the \( B \)-solution, relative to the alteration in numerical cost resulting from this coarse-graining. We note the special role of the value \( \gamma = 4 = \sqrt{D_A / D_B} \) which can be considered a good starting point for a practical implementation of the mgRDME as it balances the jump rates of the two species and enhances the model's efficiency.

Remark 2 We note that the scenario would change if we were to define an error metric that incorporates the distance within the A-population. This is evident because  \( \|\mathbf{A}-\bar{\mathbf{A}}(K_{B},\gamma)\|_{1} \)  rises with  \( \gamma \) , irrespective of  \( K_{B} \) .

In the examined first-order reaction system, only the production of A depends on location, whereas the other reactions occur independently of the molecules' spatial position. A stronger effect of the grid size on the reaction dynamics is expected for a second-order reaction system, which is studied next.

## 3 Fast-Slow Morphogen Gradient with Dimerization

We consider the fast-slow morphogen gradient system of Sect. 2, replacing the conversion reaction \( A \to B \) by the second-order reaction of dimerization \( A + A \to B \), i.e., the chemical reactions (1) are changed to (2), while the diffusion part of the reaction-diffusion model remains consistent with that described in the previous section. Again,

Springer

6

Page 16 of 31

R. Erban, S. Winkelmann

![img-9.jpeg](img-9.jpeg)

Fig. 5 Error and cost depending on the method. Error \(\mathrm{err}(K_B, \gamma)\) as defined in (36) and numerical cost \(c(K_B, \gamma)\) as defined in (40) for different values of \(K_B\) and \(\gamma = K_B / K_A\) (given by labels of the form \([K_B \gamma]\)). The colored symbols are for orientation: blue pluses for the standard RDME, red stars for the generalized mgRDME with \(K_B = 120\) and green circles for the generalized mgRDME with \(K_B = 60\). Parameters are given in (18) and \(L = 2\) (Color figure online)

we choose the one-dimensional domain \([0, L]\). The associated modeling approaches are introduced in Sect. 3.1 and compared in Sect. 3.2.

### 3.1 Modeling Approaches for the Dimerization System

Considering the compartment-based models, the diffusion operators are the same as in Sect. 2, i.e., for the standard RDME model it is given by \(\mathcal{D}\) as defined in (11), and for the generalized mgRDME model it is given by \(\tilde{\mathcal{D}}\) as defined in (21). In contrast, the reaction operators now explicitly depend on the grid size \(h\) via the second-order reaction \(A + A \to B\). For the standard compartment-based model, given by the RDME \(\frac{\partial p}{\partial t} (\mathbf{n}, \mathbf{m}, t) = (\mathcal{D} + \mathcal{R}) p(\mathbf{n}, \mathbf{m}, t)\), it reads

\[
\begin{array}{l} \mathcal {R} f (\mathbf {n}, \mathbf {m}) := k _ {1} \left\{f \left(\mathcal {O} _ {1} ^ {-} \mathbf {n}, \mathbf {m}\right) - f (\mathbf {n}, \mathbf {m}) \right\} \\ + \frac {k _ {2}}{h} \sum_ {i = 1} ^ {K} \left\{(n _ {i} + 2) (n _ {i} + 1) f (\mathcal {O} _ {i} ^ {+} \mathcal {O} _ {i} ^ {+} \mathbf {n}, \mathcal {O} _ {i} ^ {-} \mathbf {m}) - \max \{n _ {i} (n _ {i} - 1), 0 \} f (\mathbf {n}, \mathbf {m}) \right\} \\ + k _ {3} \sum_ {i = 1} ^ {K} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\} \\ \end{array}
\]

for  \( h = h_{A} = h_{B} \) . The rate constant  \( k_{2} \)  now has physical units [length]/[time], i.e., units of the ‘one-dimensional volume’ divided by time, because  \( k_{2} \)  is the rate constant of a bimolecular reaction, in contrast to  \( k_{2} \)  in (12), which has physical units [time] \( ^{-1} \) . Dividing by h in the second line stems from the standard scaling of a second-order reaction rate by the volume of the domain in which it takes place, which is given by h for the one-dimensional case under consideration. Similarly, for the generalized mgRDME model, we have

\[
\tilde {\mathcal {R}} f (\mathbf {n}, \mathbf {m}) := k _ {1} \left\{f \left(\mathcal {O} _ {1} ^ {-} \mathbf {n}, \mathbf {m}\right) - f (\mathbf {n}, \mathbf {m}) \right\}
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 17 of 31

6

\[
\begin{array}{l} + \frac {k _ {2}}{h _ {A} \gamma} \sum_ {j = 1} ^ {K _ {A}} \sum_ {i \in \mathcal {I} (j)} \left\{(n _ {j} + 2) (n _ {j} + 1) f (\mathcal {O} _ {j} ^ {+} \mathcal {O} _ {j} ^ {+} \mathbf {n}, \mathcal {O} _ {i} ^ {-} \mathbf {m}) - \max \{n _ {j} (n _ {j} - 1), 0 \} f (\mathbf {n}, \mathbf {m}) \right\} \\ + k _ {3} \sum_ {i = 1} ^ {K _ {B}} \left\{(m _ {i} + 1) f (\mathbf {n}, \mathcal {O} _ {i} ^ {+} \mathbf {m}) - m _ {i} f (\mathbf {n}, \mathbf {m}) \right\}. \\ \end{array}
\]

The factor  \( 1/\gamma \)  in the second summand results from splitting the reaction  \( A + A \rightarrow B \)  occurring in the j-th compartment into  \( |\mathcal{I}(j)| = \gamma \)  possible ones, depending on the placement of the product B:

\[
A _ {j} + A _ {j} \xrightarrow {k _ {2} / \gamma} B _ {i}, \quad j = 1, 2, \dots , K _ {A}, \quad i \in \mathcal {I} (j) := \{(j - 1) \gamma + 1, \dots , j \gamma \}, \tag {41}
\]

which conforms with the first-order case, see (20) and (22).

#### 'Ground truth' particle-based model: Brownian dynamics

As the 'ground truth' we choose particle-based dynamics given by the Doi or \(\lambda-\varrho\) model (Doi 1976a,b; Erban and Chapman 2009; Lipkova et al. 2011). Particles move in space by Brownian motion, and two particles of species \(A\) undergo dimerization \(A+A\to B\) at rate \(\lambda>0\) whenever being within a separation \(\varrho>0\), called the reaction radius.

#### Relation between rate constants for bimolecular reactions

Let \( h \) denote the grid size of the compartment-based model (given by the standard RDME), and let \( k_{2} / h^{d} \) (where \( d \) is the dimension of physical space) denote the rate for two molecules of \( A \) to react when located in the same compartment. The (standard) RDME may be seen as a formal approximation of the Doi model (Isaacson 2008), but it loses second-order reactions as \( h \to 0 \) (Isaacson and Peskin 2006; Isaacson 2009; Erban and Chapman 2009). In Isaacson (2013) a convergent RDME has been developed, which allows second-order reactions of molecules located in different compartments and thereby does not lose bimolecular reactions as \( h \to 0 \). The standard RDME may be interpreted as an asymptotic approximation of the convergent RDME for \( \varrho / h \ll 1 \) (Isaacson 2013), where \( \varrho \) is the reaction radius of the Doi model. In this case of a comparatively large compartment size, particles can be expected to only react when being located in the same box, and an adequate choice of the binding rate constant \( k_{2} \) is given by \( k_{2} = \lambda |B_{\varrho}| \), where \( |B_{\varrho}| \) denotes the volume of the \( d \)-dimensional sphere of radius \( \varrho \). For our one-dimensional domain \( [0, L] \) (i.e., \( d = 1 \)) we thus choose \( k_{2} = 2\varrho \lambda \) for the reaction-rate constant of the compartment-based model with sufficiently large compartment size \( h \). Note that, on the other hand, \( h \) has to be small, \( h \ll L \), to ensure an appropriate level of spatial resolution (Erban and Chapman 2020).

### 3.2 Model Comparison

In contrast to the first-order reaction-diffusion system of Sect. 2, analytical insights for the dimerization system are relatively limited, so our studies solely rely on computational simulations. Figure 6 shows the steady-state gradients of molecules of A and B, estimated from long-term stochastic simulations (using the Gillespie stochastic

Springer

6

Page 18 of 31

R. Erban, S. Winkelmann

![img-10.jpeg](img-10.jpeg)

Fig. 6 Steady-state solution for the morphogen system with the dimerization reaction. Comparison between Brownian dynamics (solid lines) and generalized mgRDME model (dots/bars) for dimerization system (2). The parameters are given in (42), together with L = 1,  \( K_{A} = 15 \)  and  \( K_{B} = 60 \)  (Color figure online)

simulation algorithm for the generalized mgRDME and temporal discretization for the Brownian dynamics model). The parameters are chosen as

\[
k _ {1} = 5 0, \lambda = 5, \varrho = 0. 0 2, k _ {2} = 2 \lambda \varrho = 0. 2, k _ {3} = 2, D _ {A} = 0. 1 6, D _ {B} = 0. 0 1. \tag {42}
\]

We observe a close agreement of the generalized compartment-based model with the particle-based dynamics for  \( K_{A} = 15 \)  and  \( K_{B} = 60 \) . In contrast, for the standard RDME the agreement is worse, as we will see in the next section.

#### 3.2.1 Multi-Objective Optimization

For the spatial error in the distribution of molecules of B, we define the function  \( \text{err}(K_{B}, \gamma) \)  in analogy to (36), that is, we choose the  \( L^{1} \) -distance between the steady-state solutions of the compartment-based model and of the particle-based model. Let  \( \tilde{A}_{total} \)  and  \( \tilde{B}_{total} \)  be the total number of molecules of A and B, respectively, of the compartment-based dimerization system at steady state. In contrast to the first-order reaction system of Sect. 2, we lack analytical expression for these quantities and can only estimate them via long-term numerical simulation. The cost function is defined similarly to (40), using the total propensity for diffusion and reaction events at steady state:

\[
c (K _ {B}, \gamma) := 2 \frac {D _ {A}}{h _ {A} ^ {2}} \tilde {A} _ {\text { total }} + 2 \frac {D _ {B}}{h _ {B} ^ {2}} \tilde {B} _ {\text { total }} + k _ {1} + k _ {2} \tilde {A} _ {\text { total }} (\tilde {A} _ {\text { total }} - 1) + k _ {3} \tilde {B} _ {\text { total }}, \tag {43}
\]

which is quadratic in the total number of molecules of A at steady state. Figure 7 illustrates the error and cost values across different combinations of  \( K_{B} \)  and  \( \gamma \) . We see that the error is non-monotone in  \( \gamma \) , which differs from the behavior seen in Fig. 5 for the first-order reaction-diffusion system. It seems that there is an optimal  \( K_{A} \approx 15 = 120/8 = 60/4 \) , where the error in the spatial distribution of molecules of B is minimal. Next, we study the relationship between the parameter  \( K_{A} \)  and the error in the A-population.

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 19 of 31

6

![img-11.jpeg](img-11.jpeg)

Fig. 7 Error and cost for dimerization. Cost \( c(K_B, \gamma) \) as defined in (43) and error \( \text{err}(K_B, \gamma) \) given by (36) for different values of \( K_B \) and \( \gamma = K_B / K_A \) (given by labels of the form \( [K_B \gamma] \)). The parameters are given in (42) and \( L = 1 \). The results are calculated as averages over long-term simulations of length \( T = 5000 \) (Color figure online)

#### 3.2.2 Error Analysis for Species A

Given the steady-state solution \(\bar{\mathbf{A}} = (\bar{A}_1, \bar{A}_2, \dots, \bar{A}_{K_A})\) of the generalized mgRDME model for the dimerization system, we define

\[
\bar {a} _ {h _ {A}} (x) := \bar {A} _ {i} / h _ {A} \quad \text { for } x \in [ (i - 1) h _ {A}, i h _ {A}), \quad i = 1, 2, \dots , K _ {A}, \tag {44}
\]

as an approximation of the position-dependent steady-state concentration. Let further  \( \bar{a}(x) \)  be the steady-state solution of the particle-based Brownian-dynamics dimerization model (these functions can be approximately determined via numerical simulation). We define the spatial error in A in analogy to (36):

\[
\operatorname{err} _ {A} (K _ {A}) := | | \bar {a} - \bar {a} _ {h} | | _ {L ^ {1}} = \int_ {0} ^ {L} | \bar {a} (x) - \bar {a} _ {h} (x) | \mathrm{d} x. \tag {45}
\]

Moreover, we consider the difference in the total number of molecules of \( A \) at steady state:

\[
\operatorname{err} _ {\text { total }} (K _ {A}) := \left| \bar {A} _ {\text { total }} - \int_ {0} ^ {L} \bar {a} (x)   \mathrm{d} x \right|. \tag {46}
\]

As opposed to the first-order system, this quantity can here be non-zero. Figure 8 shows a non-monotone relationship between the number \( K_{A} \) of compartments and the error \( \mathrm{err}_{\mathrm{total}}(K_A) \) in the total number of molecules of \( A \). The non-monotony is due to a too small grid size \( h_A \) for large \( K_{A} \) which contrasts with the condition \( \varrho / h \ll 1 \) necessary for a good approximation of second-order reactions, see also the paragraph on the relation between the rates on page 17. The spatial error \( \mathrm{err}_A(K_A) \) in \( A \), however, decreases monotonically with increasing \( K_{A} \) because a higher level of spatial resolution dominates the effect of losing second-order reactions for the parameter values under consideration. The error in \( B \) (also depicted in Fig. 8b) results from a combination of the two errors in \( A \) and minimizes for \( K_{A} \approx 16 \), which is in consistency with the observations from Fig. 7. The standard compartment-based model would mean to choose \( K_{A} = K_{B} = 60 \), inducing a comparatively large error both in

Springer

6

Page 20 of 31

R. Erban, S. Winkelmann

![img-12.jpeg](img-12.jpeg)

(a) Error in \(A\) for \(K_{A}\in \{2,4,\dots ,20\}\)

![img-13.jpeg](img-13.jpeg)

(b) Error in \(A\) and \(B\) for \(K_{A} \in \{10, 12, \dots, 20\}\)

Fig. 8 Error in dependence on \( K_{A} \). a Spatial error \( \mathrm{err}_A(K_A) \) given by (45) (blue) and difference \( \mathrm{err}_{\mathrm{total}}(K_A) \) in the total number of molecules of \( A \) defined in (46) (red) for different values of \( K_{A} \). b zoom-in, additionally contains the spatial error in \( B \), \( \mathrm{err}(K_B, \gamma) \) defined in (36), for fixed \( K_{B} = 60 \) and \( \gamma = K_{B} / K_{A} \). Parameters given in (42), \( L = 1 \). Long-term simulation of length \( T = 10^{5} \) (Color figure online)

the steady-state distribution of molecules of \( B \) and in the total number of molecules of \( A \).

We conclude that the generalized mgRDME model clearly outperforms the standard RDME model for the fast-slow morphogen gradient system with dimerization. For practical implementation, the preliminary choice of \(\gamma = \sqrt{D_A / D_B}\) must be carefully adapted to ensure the adequate reproduction of the second-order reaction \(A + A\rightarrow B\). Specifically, the compartment size \(h_A\) should be chosen based on the reaction radius \(\varrho\) and the rate constant \(\lambda\) of the particle-based Brownian dynamics model, in conjunction with the diffusion constant \(D_{A}\), but independently of the spatial resolution for \(B\). This stands in clear contrast to the first-order chemical reaction system discussed in Sect. 2, where it is more appropriate to adjust the compartment sizes of the two species relative to each other.

Remark 3 (Second-order moments for dimerization) Unlike the first-order reaction-diffusion model discussed in Sect. 2, the variance of the particle numbers does here not align with the mean. Consequently, one could also take the variance as a value for comparing the models and defining an error. A numerical comparison of mean and variance for the different model types in application to the dimerization system is given in Appendix 5.3.

## 4 Discussion

We examined the mgRDME as a generalized compartment-based stochastic model for reaction-diffusion kinetics, allowing adapted compartment sizes for different chemical species involved in the system (Cao and Erban 2014). When employed in the context of morphogen gradient formation, the mgRDME is useful because the signaling molecules A typically diffuse at a significantly higher rate than the morphogen molecules B, and this diffusion rate is a deciding factor for the choice of compartment sizes. In comparison to the 'ground truth' particle-based Brownian dynamics simulations, the mgRDME outperforms the standard RDME in terms of accuracy and

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 21 of 31

6

reduction in numerical cost, both for the first-order reaction network (1) and for the considered dimerization system (2).

Interestingly, however, the structural dependence of error and cost on the compartment sizes is different for the two settings. In Sect. 2, the optimal ratio between grid sizes hinges on the relationship between the diffusion rates. Conversely, in the system involving second-order dimerization in Sect. 3, an optimal grid size for species \( A \) can be identified (in dependence on its diffusion constant) regardless of the grid size allocated to \( B \). In particular, the approximation error may be reduced for the dimerization system by choosing a smaller compartment size for \( B \) while maintaining the fixed optimal compartment size for \( A \) undergoing dimerization. This goes in line with the observations for standard RDMEs, where very small grid sizes result in the elimination of the second-order reactions in three-dimensional simulations. The mgRDME circumvents this issue by allowing larger compartment sizes for species involved in second-order reactions. At the same time, the increase in numerical cost caused by the higher spatial resolution is limited, since it is selectively applied solely to the slow species.

Our work shows that the mgRDME is a valuable framework for studying reaction-diffusion systems with multiple scales in diffusion speed of involved molecules. The multi-grid approach enables the mgRDME to more effectively accommodate the inherent conditions of the system under consideration, surpassing the adaptability of the standard RDME. The models with different diffusion constants naturally take into account variations in the sizes of simulated biomolecules, which can range from small calcium ions to relatively large vesicles in neurotransmission dynamics (Ernst et al. 2022, 2023), or from relatively small \(G\)-actin monomers to larger \(F\)-actin filaments in simulations of actin dynamics (Zhuravlev et al. 2010; Dobramysl et al. 2016b). Other applications of mgRDME include pattern formation based on Turing instability which requires different diffusion constants (Turing 1952).

In this work, we have solely considered dynamics in one-dimensional domains – a restriction which is well justified in the context of applications to developmental biology, where a modeller is interested in the pattern formation along one body axis, see, for example, dorsal-ventral patterning studies of the Drosophila embryo (Shimmi and O'Connor 2003; Shimmi et al. 2005). However, the ideas behind mgRDME can be extended to higher-dimensional models, where unstructured meshes are used to discretize complex geometries of biological problems to better resolve curved inner and outer boundaries than by Cartesian meshes (Engblom et al. 2009; Isaacson and Zhang 2018). The mgRDME framework can use different unstructured meshes to describe different chemical species. In two-dimensional domains, voxels are typically triangles or squares, while in three-dimensional domains, they are tetrahedra or cubes (Engblom et al. 2009; Earnest et al. 2015). The reaction rates in mgRDME are constructed based on the shared area or volume between voxels from different grids, as described in Remark 1, where we note that voxels may partially overlap, see Eq. (25). The relative size of different meshes (parameterized by \(\gamma\) in this paper) can be chosen to optimize the number of spatial diffusive jumps for species with varying diffusion constants. Additionally, different species can be modeled at different resolutions to improve simulation accuracy, with voxel sizes optimized to enhance the description of bimolecular reactions. In applications where high accuracy is only needed in a subset

Springer

6

Page 22 of 31

R. Erban, S. Winkelmann

of the simulation domain (Dobramysl et al. 2016a, b), unstructured meshes also offer the advantage to refine the mesh in the corresponding subset of the computational domain. In particular, to model complex biological systems, the mgRDME approach (of using different meshes for different chemical species) can be combined with using different mesh resolutions (for the same chemical species) in different parts of the computational domain (Kang and Erban 2019).

While our analysis has been based on comparison of steady-state morphogen gradients, another way to compare the RDME and mgRDME could be by comparing the first collision and mean reaction times in these models (Li et al. 2018), which can be achieved by analyzing the random walks on lattices (Montroll 1969). Also a comparison of time-dependent solutions (in contrast to steady-state distributions investigated here) is another interesting topic for future research on mgRDME models of morphogen gradients (Saunders and Howard 2009). Some approaches to compare the dynamics prior to steady state are discussed in Appendix 5.4.

## 5 Appendix

### 5.1 Time-Dependent Solution of the Brownian Dynamics Model

Let us consider that initially there are no molecules of \( A \) and \( B \) in the system, i.e.

\[
a (x, 0) \equiv b (x, 0) \equiv 0. \tag {47}
\]

To solve the PDE (27) with initial condition (47) and boundary conditions (29), we apply the Fourier transform

\[
\widehat {a} (\xi , t) = \int_ {- \infty} ^ {\infty} a (x, t) e ^ {- i \xi x} d x.
\]

We get

\[
\frac {\partial \widehat {a}}{\partial t} (\xi , t) = - (D _ {A} \xi^ {2} + k _ {2}) \widehat {a} (\xi , t) + 2 k _ {1}.
\]

Solving this ODE with the initial condition  \( \widehat{a}(\xi,0)=0 \) , we obtain

\[
\widehat {a} (\xi , t) = \frac {2 k _ {1}}{D _ {A} \xi^ {2} + k _ {2}} - \frac {2 k _ {1} \exp \left[ - \left(D _ {A} \xi^ {2} + k _ {2}\right) t \right]}{D _ {A} \xi^ {2} + k _ {2}}. \tag {48}
\]

Since the Fourier transform of the convolution is the product of Fourier transforms, we can calculate the Fourier inverse of (48) as follows

\[
a (x, t) = \bar {a} (x) - \frac {\exp [ - k _ {2} t ]}{2 \sqrt {D _ {A} t \pi}} \int_ {- \infty} ^ {\infty} \bar {a} (y) \exp \left[ - \frac {(x - y) ^ {2}}{4 D _ {A} t} \right] d y,
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 23 of 31

6

where \(\tilde{a}(x)\) is given by (33). In total we get

\[
\begin{array}{l} a (x, t) = \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | x | \right] \\ - \frac {\exp [ - k _ {2} t ]}{2 \sqrt {D _ {A} t \pi}} \int_ {- \infty} ^ {\infty} \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | y | \right] \exp \left[ - \frac {(x - y) ^ {2}}{4 D _ {A} t} \right] d y \tag {49} \\ \end{array}
\]

for the average number of \(A\)-molecules in the interval \([x, x + \mathrm{d}x]\) at time \(t\).

### 5.2 'Ground Truth' Solution: The Effect of the Boundary

Here, we derive the 'ground truth' solution of Sect. 2.3 for zero-flux boundary conditions. To do that, we will solve the ODEs (30-31) on the interval \([-L, L]\) with zero-flux boundary conditions

\[
\frac {\mathrm{d} \tilde {a}}{\mathrm{d} x} (- L) = \frac {\mathrm{d} \tilde {a}}{\mathrm{d} x} (L) = \frac {\mathrm{d} \tilde {b}}{\mathrm{d} x} (- L) = \frac {\mathrm{d} \tilde {b}}{\mathrm{d} x} (L) = 0. \tag {50}
\]

Then the solution to the steady-state Eq. (30) is

\[
\tilde {a} (x) = \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} | x | \right] \frac {1 + \exp \left[ 2 \sqrt {\frac {k _ {2}}{D _ {A}}} (| x | - L) \right]}{1 - \exp [ - 2 L \sqrt {k _ {2} / D _ {A}} ]},
\]

which for \(L \to \infty\) agrees with (33), independently of the parameter values.

### 5.3 Second-Order Moments for Dimerization

Means and variances of the numbers of molecules of A and B for the dimerization system are plotted in Fig. 9, sampled from long-term simulations of the (standard) compartment-based model and of the particle-based Brownian dynamics model, choosing boxes of size h = 1/15. We observe a good agreement of the two model types.

### 5.4 Dynamics Prior to Steady State

To study the dynamics prior to steady state we consider the average positions of the particles in dependence on time. These investigations are done for the first-order morphogen gradient system introduced in Sect. 2.

Springer

6

Page 24 of 31

R. Erban, S. Winkelmann

![img-14.jpeg](img-14.jpeg)

(a) mean \(A\)

![img-15.jpeg](img-15.jpeg)

(b) mean \(B\)

![img-16.jpeg](img-16.jpeg)

(c) variance \(A\)

![img-17.jpeg](img-17.jpeg)

(d) variance \(B\)

Fig. 9 Comparison of mean and variance for dimerization. Mean/variance of the number of molecules of a/c species A and b/d species B at steady state estimated by long-term simulation (T = 2000) of the compartment-based model (red circles) and of the particle-based Brownian dynamics model (grey bars). Rate constants given in (42), L = 1,  \( K_{A} = K_{B} = 15 \)  ( \( \gamma = 1 \) , h = 1/15) (Color figure online)

#### 5.4.1 Average Positions for 'Ground Truth' Solution

We define the moments

\[
\langle 1 \rangle_ {a} (t) := \int_ {0} ^ {\infty} a (x, t) \mathrm{d} x, \quad \langle x \rangle_ {a} (t) := \int_ {0} ^ {\infty} x a (x, t) \mathrm{d} x, \tag {51}
\]

\[
\langle 1 \rangle_ {b} (t) := \int_ {0} ^ {\infty} b (x, t) \mathrm{d} x, \quad \langle x \rangle_ {b} (t) := \int_ {0} ^ {\infty} x b (x, t) \mathrm{d} x. \tag {52}
\]

Multiplying Eqs. (27-28) by 1 and \( x \), respectively, and integrating over \( x \in [0, \infty) \), we obtain the ODEs

\[
\frac {\mathrm{d} \langle 1 \rangle_ {a}}{\mathrm{d} t} = - k _ {2} \langle 1 \rangle_ {a} + k _ {1}, \tag {53}
\]

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 25 of 31

6

\[
\frac {\mathrm{d} \langle x \rangle_ {a}}{\mathrm{d} t} = D _ {A} a (0, t) - k _ {2} \langle x \rangle_ {a}, \tag {54}
\]

\[
\frac {\mathrm{d} \langle 1 \rangle_ {b}}{\mathrm{d} t} = k _ {2} \langle 1 \rangle_ {a} - k _ {3} \langle 1 \rangle_ {b}, \tag {55}
\]

\[
\frac {\mathrm{d} \langle x \rangle_ {b}}{\mathrm{d} t} = D _ {B} b (0, t) + k _ {2} \langle x \rangle_ {a} - k _ {3} \langle x \rangle_ {b}. \tag {56}
\]

Substituting \( x = 0 \) into Eq. (49), we have

\[
\begin{array}{l} a (0, t) = \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} - \frac {\exp [ - k _ {2} t ]}{\sqrt {D _ {A} t \pi}} \int_ {0} ^ {\infty} \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \exp \left[ - \sqrt {\frac {k _ {2}}{D _ {A}}} y - \frac {y ^ {2}}{4 D _ {A} t} \right] d y \\ = \frac {k _ {1}}{\sqrt {D _ {A} k _ {2}}} \operatorname{erf} \left[ \sqrt {k _ {2} t} \right]. \tag {57} \\ \end{array}
\]

Substituting (57) into Eq. (54), we get

\[
\frac {\mathrm{d} \langle x \rangle_ {a}}{\mathrm{d} t} = \frac {k _ {1} \sqrt {D _ {A}}}{\sqrt {k _ {2}}} \operatorname{erf} \left[ \sqrt {k _ {2} t} \right] - k _ {2} \langle x \rangle_ {a}.
\]

Using the initial condition \(\langle x\rangle_a(0) = 0\), we obtain

\[
\begin{array}{l} \langle x \rangle_ {a} (t) = \exp [ - k _ {2} t ] \int_ {0} ^ {t} \exp [ k _ {2} \tau ] \frac {k _ {1} \sqrt {D _ {A}}}{\sqrt {k _ {2}}} \operatorname{erf} [ \sqrt {k _ {2} \tau} ] d \tau \\ = \frac {k _ {1} \sqrt {D _ {A}}}{k _ {2} \sqrt {k _ {2}}} \operatorname{erf} \left[ \sqrt {k _ {2} t} \right] - \exp \left[ - k _ {2} t \right] \frac {2 k _ {1} \sqrt {D _ {A} t}}{k _ {2} \sqrt {\pi}}. \tag {58} \\ \end{array}
\]

Using the initial condition \(\langle 1\rangle_{a}(0) = 0\) in Eq. (53), we obtain

\[
\langle 1 \rangle_ {a} (t) = \frac {k _ {1}}{k _ {2}} \left(1 - \exp [ - k _ {2} t ]\right) \tag {59}
\]

for the average total number of \( A \)-particles at time \( t \). Consequently, the average position \( x_{a}(t) \) of \( A \)-particles as a function of time \( t > 0 \) is

\[
x _ {a} (t) := \frac {\langle x \rangle_ {a} (t)}{\langle 1 \rangle_ {a} (t)} = \frac {\sqrt {D _ {A} \pi} \operatorname{erf} [ \sqrt {k _ {2} t} ] - \exp [ - k _ {2} t ] 2 \sqrt {D _ {A} k _ {2} t}}{(1 - \exp [ - k _ {2} t ]) \sqrt {k _ {2} \pi}}. \tag {60}
\]

We note that this ratio is independent of our choice of the rate constant \( k_{1} \).

Substituting (59) into Eq. (55), we get

\[
\frac {\mathrm{d} \langle 1 \rangle_ {b}}{\mathrm{d} t} = k _ {1} (1 - \exp [ - k _ {2} t ]) - k _ {3} \langle 1 \rangle_ {b}.
\]

Springer

6

Page 26 of 31

R. Erban, S. Winkelmann

![img-18.jpeg](img-18.jpeg)

(a) Average number of particles

![img-19.jpeg](img-19.jpeg)

(b) Average positions

Fig. 10 'Ground truth' dynamics prior to steady state. a Average total number of particles \(\langle 1\rangle_{a}(t)\) and \(\langle 1\rangle_{b}(t)\) given by (59) and (61), respectively. b Average positions \(x_{a}(t)\) and \(x_{b}(t)\) depending on time given by (60) and (62). Rate constants given in (18) (Color figure online)

Using the initial condition \(\langle 1\rangle_b(0) = 0\), we obtain

\[
\begin{array}{l} \langle 1 \rangle_ {b} (t) = \exp [ - k _ {3} t ] \int_ {0} ^ {t} \exp [ k _ {3} \tau ] k _ {1} (1 - \exp [ - k _ {2} \tau ]) d \tau \\ = \frac {k _ {1}}{k _ {3}} + \frac {k _ {1} k _ {2} \exp [ - k _ {3} t ]}{(k _ {3} - k _ {2}) k _ {3}} + \frac {k _ {1} \exp [ - k _ {2} t ]}{k _ {2} - k _ {3}} \tag {61} \\ \end{array}
\]

for the average total number of \( B \)-particles at time \( t \). Substituting (58) into Eq. (56), we get

\[
\frac {\mathrm{d} \langle x \rangle_ {b}}{\mathrm{d} t} = D _ {B} b (0, t) + \frac {k _ {1} \sqrt {D _ {A}}}{\sqrt {k _ {2}}} \operatorname{erf} \left[ \sqrt {k _ {2} t} \right] - \exp [ - k _ {2} t ] \frac {2 k _ {1} \sqrt {D _ {A} t}}{\sqrt {\pi}} - k _ {3} \langle x \rangle_ {b}.
\]

The average position  \( x_{b}(t) \)  of B-particles at time t is given by

\[
x _ {b} (t) = \frac {\langle x \rangle_ {b} (t)}{\langle 1 \rangle_ {b} (t)}. \tag {62}
\]

Fig. 10 shows the average total number of particles and their average positions depending on time for the ‘ground truth’ model.

#### 5.4.2 Time-Dependent Solution for the Compartment-Based Model

In analogy to the densities \( a(x, t) \) and \( b(x, t) \) defined for the particle-based model, let \( a_i(t) \) and \( b_i(t) \) denote the average number of particles of species \( A \) and \( B \) in box \( i \) at time \( t \), respectively. We set \( \mathbf{a}(t) = (a_i(t))_{i=1,\dots,K_A} \) and \( \mathbf{b}(t) = (b_i(t))_{i=1,\dots,K_B} \) and

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 27 of 31

6

![img-20.jpeg](img-20.jpeg)

(a) Difference depending on \(\gamma\) \((K_B = 120)\)

![img-21.jpeg](img-21.jpeg)

(b) Difference depending on \(K_{B}\) \((\gamma = 1)\)

Fig. 11 Time-dependent differences of average positions. a Differences between average positions \( x_B^{(\gamma)}(t) \) of mgRDME, see (69), and average positions \( x_b(t) \) of the 'ground truth' solution, see (62), for various \( \gamma \) and fixed \( K_B = 120 \). b Differences between average positions of the RDME (fixed \( \gamma = 1 \)) and the 'ground truth' solution for various \( K_A = K_B \). Rate constants given in (18) (Color figure online)

have

\[
\frac {\mathrm{d} \mathbf {a}}{\mathrm{d} t} = \left(\frac {D _ {A}}{h _ {A} ^ {2}} S _ {A} - k _ {2} I _ {A}\right) \mathbf {a} + k _ {1} \mathbf {e} _ {1}, \tag {63}
\]

\[
\frac {\mathrm{d} \mathbf {b}}{\mathrm{d} t} = \left(\frac {D _ {B}}{h _ {B} ^ {2}} S _ {B} - k _ {3} I _ {B}\right) \mathbf {b} + \frac {k _ {2}}{\gamma} M \mathbf {a}, \tag {64}
\]

for  \( I_{A}, S_{A} \in R^{K_{A} \times K_{A}}, I_{B}, S_{B} \in R^{K_{B} \times K_{B}}, e_{1} \in R^{K_{A}} \)  and  \( M \in R^{K_{B} \times K_{A}} \)  defined in Sect. 2.2. We note that the standard RDME is included by setting  \( K_{A} = K_{B} \) . The ODEs (63)-(64) are the analogue to Eqs. (27–28) of the ‘ground truth’ dynamics. Together with  \( \mathbf{a}(0) = 0 \)  we obtain

\[
\mathbf {a} (t) = W _ {A} ^ {- 1} \left[ \exp \left(W _ {A} t\right) - I _ {A} \right] k _ {1} \mathbf {e} _ {1} \tag {65}
\]

for

\[
W _ {A} := \frac {D _ {A}}{h _ {A} ^ {2}} S _ {A} - k _ {2} I _ {A} \in \mathbb {R} ^ {K _ {A} \times K _ {A}}.
\]

Moreover, using \(\mathbf{b}(0) = 0\) and defining \(W_{B} := (D_{B} / h_{B}^{2})S_{B} - k_{3}I_{B} \in \mathbb{R}^{K_{B} \times K_{B}}\), we obtain (Perko 2013)

\[
\begin{array}{l} \mathbf {b} (t) = \frac {k _ {2}}{\gamma} \int_ {0} ^ {t} \exp (W _ {B} (t - s)) M \mathbf {a} (s) d s \\ \stackrel {(6 5)} {=} \frac {k _ {2}}{\gamma} \int_ {0} ^ {t} \exp (W _ {B} (t - s)) M W _ {A} ^ {- 1} \left[ \exp \left(W _ {A} s\right) - I _ {A} \right] k _ {1} \mathbf {e} _ {1} d s \tag {66} \\ \end{array}
\]

Springer

6

Page 28 of 31

R. Erban, S. Winkelmann

for the vector of time-dependent average numbers of \(B\)-particles in the boxes at time \(t\). We have \(\lim_{t\to \infty}\mathbf{a}(t) = \overline{\mathbf{A}}\) and \(\lim_{t\to \infty}\mathbf{b}(t) = \overline{\mathbf{B}}\) for the steady-state solutions \(\overline{\mathbf{A}},\overline{\mathbf{B}}\) given in (16) and (24).

#### 5.4.3 Average Positions for the Compartment-Based Model

In analogy to the moments defined in (51)-(52), we set

\[
\langle 1 \rangle_ {A} (t) := \sum_ {i = 1} ^ {K _ {A}} a _ {i} (t), \quad \langle x \rangle_ {A} (t) := \sum_ {i = 1} ^ {K _ {A}} (i - 0. 5) h _ {A} a _ {i} (t), \tag {67}
\]

\[
\langle 1 \rangle_ {B} (t) := \sum_ {i = 1} ^ {K _ {B}} b _ {i} (t), \quad \langle x \rangle_ {B} (t) := \sum_ {i = 1} ^ {K _ {B}} (i - 0. 5) h _ {B} b _ {i} (t). \tag {68}
\]

In the definitions of  \( \langle x\rangle_{A}(t) \)  and  \( \langle x\rangle_{B}(t) \) , we choose the centers  \( (i-0.5)h_{A} \)  and  \( (i-0.5)h_{B} \)  of the boxes for the spatial value, which is justified by the assumption that the particles in each box are spatially well mixed so that their positions follow a uniform distribution.

By summing Eqs. (63–64) over the discrete index i we obtain the same ODEs for the moments  \( \langle1\rangle_{A}(t) \)  and  \( \langle1\rangle_{B}(t) \)  as in the ‘ground truth’ model; see Eqs. (53) and (55). This means that the average total number of particles  \( \langle1\rangle_{A}(t) \)  and  \( \langle1\rangle_{B}(t) \)  is again given by Eqs. (59) and (61). We calculate the moments  \( \langle x\rangle_{A}(t) \)  and  \( \langle x\rangle_{B}(t) \)  according to their definition in (67)-(68), using the solutions  \( \mathbf{a}(t) \)  and  \( \mathbf{b}(t) \)  given in (65) and (66). Given  \( \gamma = K_{B}/K_{A} \) , the mean positions  \( x_{A}^{(\gamma)}(t) \)  and  \( x_{B}^{(\gamma)}(t) \)  of A and B molecules at time t are

\[
x _ {A} ^ {(\gamma)} (t) = \frac {\langle x \rangle_ {A} (t)}{\langle 1 \rangle_ {A} (t)}, \quad x _ {B} ^ {(\gamma)} (t) = \frac {\langle x \rangle_ {B} (t)}{\langle 1 \rangle_ {B} (t)}, \tag {69}
\]

in analogy to (60) and (62), respectively.

Figure 11 shows the time-dependent differences between the average positions \( x_B^{(\gamma)}(t) \) of the mgRDME and the average positions \( x_b(t) \) of the 'ground truth' solution for various \( \gamma \) and \( K_B \). We observe that for every time \( t > 0 \), the error increases with the ratio \( \gamma = K_B / K_A \) for fixed \( K_B = 120 \), see Fig. 11a. Similarly, the error increases with decreasing values of \( K_B \), see Fig. 11b. In contrast to the results from Sect. 2.4, the approximation quality is here more robust against changes in \( K_B \); this can be explained by the fact that the reference value is given by the average positions and not by the spatial distribution of molecules of \( B \).

The mean position  \( x_{B}^{(\gamma)}(t) \)  of the compartment-based model always exceeds the mean position  \( x_{b}(t) \)  of the ‘ground truth’ solution. Especially for the initial time phase (small t), the difference is large. This is due to the fact that the production of molecules of A is located at position x = 0. In the mgRDME, the newly produced particles immediately get well mixed in the first box  \( [0, h_{A}] \) , while in the ‘ground truth’ model it takes some time for the particles to spread out in space. The mixing in the mgRDME implies a shift of the average to the right, which has a larger effect especially in the beginning, where the population is small and mostly located near the

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 29 of 31

6

origin. This effect, which also propagates to B via conversion, can be interpreted as an artificial dispersion induced by the mgRDME.

Author Contributions Both authors formulated the problem, carried out the theoretical analysis of the models and wrote the paper. SW carried out illustrative numerical simulations.

Funding This work was supported by the Engineering and Physical Sciences Research Council, grant number EP/V047469/1, awarded to Radek Erban. We further acknowledge the support of Deutsche Forschungsgemeinschaft through CRC 1114 and through Germany's Excellence Strategy–The Berlin Mathematics Research Center MATH+ (EXC-2046/1, project 390685689), awarded to Stefanie Winkelmann.

Data Availability In compliance with EPSRC's open access initiative, the data in this paper is available from: https://doi.org/10.5281/zenodo.11072450.

## Declarations

Conflict of interest Authors have no conflict of interest.

Ethical Approval Not applicable.

Open Access This article is licensed under a Creative Commons Attribution 4.0 International License, which permits use, sharing, adaptation, distribution and reproduction in any medium or format, as long as you give appropriate credit to the original author(s) and the source, provide a link to the Creative Commons licence, and indicate if changes were made. The images or other third party material in this article are included in the article's Creative Commons licence, unless indicated otherwise in a credit line to the material. If material is not included in the article's Creative Commons licence and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder. To view a copy of this licence, visit http://creativecommons.org/licenses/by/4.0/.

## References

Agbanusi IC, Isaacson SA (2014) A comparison of bimolecular reaction models for stochastic reaction-diffusion systems. Bull Math Biol 76(4):922–946

Andrews SS (2012) Spatial and stochastic cellular modeling with the Smoldyn simulator. Methods and Protocols, Bacterial Molecular Networks, pp 519–542

Andrews SS, Addy NJ, Brent R, Arkin AP (2010) Detailed simulations of cell biology with Smoldyn 2.1. PLoS Comput Biol 6(3):1000705

Ben-Avraham D, Redner S (1986) Kinetics of n-species annihilation: mean-field and diffusion-controlled limits. Phys Rev A 34(1):501

Bergmann S, Sandler O, Sberro H, Shnider S, Schejter E, Shilo B, Barkai N (2007) Pre-steady-state decoding of the Bicoid morphogen gradient. PLoS Biol 5(2):46

Cao Y, Erban R (2014) Stochastic Turing patterns: analysis of compartment-based approaches. Bull Math Biol 76:3051–3069

Chatterjee A, Vlachos D, Katsoulakis M (2004) Spatially adaptive lattice coarse-grained Monte Carlo simulations for diffusion of interacting molecules. J Chem Phys 121:11420–11431

Dai J, Seider W, Sinno T (2008) Coarse-grained lattice kinetic Monte Carlo simulation of systems of strongly interacting particles. J Chem Phys 128:194705

Denizot A, Arizono M, Nägerl U, Soula H, Berry H (2019) Simulation of calcium signaling in fine astrocytic processes: effect of spatial properties on spontaneous activity. PLoS Comput Biol 15(8):1006795

Dobramysl U, Rüdiger S, Erban R (2016) Particle-based multiscale modeling of calcium puff dynamics. Multiscale Modelling Simulation 14(3):997–1016

Dobramysl U, Papoian G, Erban R (2016) Steric effects induce geometric remodeling of actin bundles in filopodia. Biophys J 110:2066–2075

Springer

6

Page 30 of 31

R. Erban, S. Winkelmann

Doi M (1976) Stochastic theory of diffusion-controlled reaction. J Phys A: Math Gen 9(9):1479

Doi M (1976) Second quantization representation for classical many-particle system. J Phys A: Math Gen 9(9):1465

Earnest TM, Lai J, Chen K, Hallock MJ, Williamson JR, Luthey-Schulten Z (2015) Toward a whole-cell model of ribosome biogenesis: kinetic modeling of SSU assembly. Biophys J 109(6):1117–1135

Engblom S, Ferm L, Hellander A, Lötstedt P (2009) Simulation of stochastic reaction-diffusion processes on unstructured meshes. SIAM J Sci Comput 31(3):1774–1797

Erban R, Chapman SJ (2007) Reactive boundary conditions for stochastic simulations of reaction-diffusion processes. Phys Biol 4(1):16–28

Erban R, Chapman SJ (2009) Stochastic modelling of reaction-diffusion processes: algorithms for bimolecular reactions. Phys Biol 6(4):046001

Erban R, Chapman SJ (2020) Stochastic Modelling of Reaction-Diffusion Processes. Cambridge University Press, Cambridge

Erban R, Flegg M, Papoian G (2014) Multiscale stochastic reaction-diffusion modelling: application to actin dynamics in filopodia. Bull Math Biol 76(4):799–818

Ernst A, Schütte C, Sigrist SJ, Winkelmann S (2022) Variance of filtered signals: characterization for linear reaction networks and application to neurotransmission dynamics. Math Biosci 343:108760

Ernst A, Unger N, Schütte C, Walter AM, Winkelmann S (2023) Rate-limiting recovery processes in neurotransmission under sustained stimulation. Math Biosci 362:109023

Fange D, Elf J (2006) Noise-induced Min phenotypes in E. coli. PLoS Comput Biol 2(6):637–648

Flegg M, Chapman J, Erban R (2012) The two-regime method for optimizing stochastic reaction-diffusion simulations. J R Soc Interface 9(70):859–868

Flegg M, Hellander S, Erban R (2015) Convergence of methods for coupling of microscopic and mesoscopic reaction-diffusion simulations. J Comput Phys 289:1–17

Gillespie DT (1977) Exact stochastic simulation of coupled chemical reactions. J Phys Chem 81(25):2340–2361

Gillespie DT (2007) Stochastic simulation of chemical kinetics. Annu Rev Phys Chem 58:35–55

Gregor T, Wieschaus E, McGregor A, Bialek W, Tank D (2007) Stability and nuclear dynamics of the bicoid morphogen gradient. Cell 130(1):141–152

Hellander S, Petzold L (2016) Reaction rates for a generalized reaction-diffusion master equation. Phys Rev E 93(1):013307

Hellander S, Hellander A, Petzold L (2012) Reaction-diffusion master equation in the microscopic limit. Phys Rev E 85(4):042901

Hellander S, Hellander A (2020) Hierarchical algorithm for the reaction-diffusion master equation. Journal of Chemical Physics 152(3)

Isaacson SA (2013) A convergent reaction-diffusion master equation. Journal of Chemical Physics 139(5)

Isaacson SA (2008) Relationship between the reaction-diffusion master equation and particle tracking models. J Phys A: Math Theor 41(6):065003

Isaacson SA (2009) The reaction-diffusion master equation as an asymptotic approximation of diffusion to a small target. SIAM J Appl Math 70(1):77–111

Isaacson SA, Peskin CS (2006) Incorporating diffusion in complex geometries into stochastic chemical kinetics simulations. SIAM J Sci Comput 28(1):47–74

Isaacson SA, Zhang Y (2018) An unstructured mesh convergent reaction-diffusion master equation for reversible reactions. J Comput Phys 374:954–983

Iyer K, Prabhakara C, Mayor S, Rao M (2023) Cellular compartmentalisation and receptor promiscuity as a strategy for accurate and robust inference of position during morphogenesis. Elife 12:79257

Jahnke T, Huisinga W (2007) Solving the chemical master equation for monomolecular reaction systems analytically. J Math Biol 54:1–26

Kang H, Erban R (2019) Multiscale stochastic reaction-diffusion algorithms combining Markov chain models with stochastic partial differential equations. Bull Math Biol 81(8):3185–3213

Kang K, Redner S (1984) Scaling approach for the kinetics of recombination processes. Phys Rev Lett 52(12):955

Kang H, Zheng L, Othmer H (2012) A new method for choosing the computational cell in stochastic reaction-diffusion systems. J Math Biol 65:1017–1099

Kicheva A, Pantazis P, Bollenbach T, Kalaidzidis Y, Bittig T, Julicher F, González-Gaitán M (2007) Kinetics of morphogen gradient formation. Science 315(5811):521–525

Springer

Multi-Grid Stochastic Reaction-Diffusion Modelling

Page 31 of 31 6

Kolomeisky AB (2011) Formation of a morphogen gradient: acceleration by degradation. *J Phys Chem Lett* 2(13):1502–1505Li F, Chen M, Erban R, Cao Y (2018) Reaction time for trimolecular reactions in compartment-based reaction-diffusion models. *J Chem Phys* 148(20):204108Li F, Cao Y (2012) Multiscale discretization for reaction diffusion systems. In: *Proceedings of the International Conference on Bioinformatics & Computational Biology (BIOCOMP)*, p. 242Lipkova J, Zygalakis K, Chapman J, Erban R (2011) Analysis of Brownian dynamics simulations of reversible bimolecular reactions. *SIAM J Appl Math* 71(3):714–730Montefusco A, Schütte C, Winkelmann S (2023) A route to the hydrodynamic limit of a reaction-diffusion master equation using gradient structures. *SIAM J Appl Math* 83(2):837–861Montroll E (1969) Random walks on lattices. III. calculation of first-passage times with application to exciton trapping on photosynthetic units. *J Math Phys* 10(4):753–765Murray J (2002) *Math Biol II: Spatial Models Biomed Appl*. Springer, New YorkPerko L (2013) Differential Equations and Dynamical Systems, vol 7. Springer, New YorkRogers K, Schier A (2011) Morphogen gradients: from generation to interpretation. *Annu Rev Cell Dev Biol* 27(1):377–407Saunders T, Howard M (2009) When it pays to rush: interpreting morphogen gradients prior to steady-state. *Phys Biol* 6(4):046020Shimmi O, O'Connor M (2003) Physical properties of Tld, Sog, Tsg and Dpp protein interactions are predicted to help create a sharp boundary in Bmp signals during dorsoventral patterning of the *Drosophila* embryo. *Development* 130(19):4673–4682Shimmi O, Umulis D, Othmer H, Connor M (2005) Facilitated transport of a Dpp/Scw heterodimer by Sog/Tsg leads to robust patterning of the *Drosophila* blastoderm embryo. *Cell* 120(6):873–886Teimouri H, Kolomeisky AB (2018) Discrete-state stochastic modeling of morphogen gradient formation. Methods and Protocols, Morphogen Gradients, pp 199–221Turing A (1952) The chemical basis of morphogenesis. *Philosophical Trans R Soc Lond* 237:37–72Winkelmann S, Schütte C (2020) Stochastic Dynamics in Computational Biology, vol 645. Springer, ChamWinkelmann S, Zonker J, Schütte C, Conrad ND (2021) Mathematical modeling of spatio-temporal population dynamics and application to epidemic spreading. *Math Biosci* 336:108619Winkelmann S, Schütte C (2016) The spatiotemporal master equation: Approximation of reaction-diffusion dynamics via Markov state modeling. *Journal of Chemical Physics* 145(21)Wolpert L, Tickle C, Arias AM (2015) Principles of Development. Oxford University Press, OxfordZhuravlev P, Der B, Papoian GA (2010) Design of active transport must be highly intricate: a possible role of myosin and ena/VASP for G-actin transport in filopodia. *Biophys J* 98:1439–1448

**Publisher's Note** Springer Nature remains neutral with regard to jurisdictional claims in published maps and institutional affiliations.

Springer
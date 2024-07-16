# CRNAS

CRNAS experiments.

This repository contains the MATLAB scripts of the experiments presented in the article "Novel Optimization Techniques for Parameter Estimation". In particular, this repository contains the following directories:
* main function: This directory contains ready-to-used CRNAS functions.
* Case study I deterministic: This directory contains codes necessary for reproduce the results presented in 'Case Study I: Cancer Drug Response Estimation - Deterministic Drug-Affected Cell Proliferation Framework'.
* Case study I stochastic: This directory contains codes necessary for reproduce the results presented in 'Case Study I: Cancer Drug Response Estimation - Stochastic Drug-Affected Cell Proliferation Framework'.
* Case study II: This directory contains codes necessary for reproduce the results presented in 'Case Study II: Heterogeneous Logistic Estimation'.

The required Matlab version is Matlab R2022a or newer.
 
 Author list: Chenyu Wu, Nuozhou Wang, Casey Garner, Kevin Leder, Shuzhong Zhang
 
 Abstract: In this paper, we introduce a new optimization algorithm that is well suited for solving parameter estimation problems. We call our new method cubic regularized Newton with affine scaling (CRNAS). In contrast to so-called first-order methods which rely solely on the gradient of the objective function, our method utilizes the Hessian of the objective. As a result it is able to focus on points satisfying the second-order optimality conditions, as opposed to first-order methods that simply converge to critical points. This is an important feature in parameter estimation problems where the objective function is often non-convex and as a result there can be many critical points making it is near impossible to identify the global minimum. An important feature of parameter estimation in mathematical models of biological systems is that the parameters are constrained by either physical constraints or prior knowledge. We use an affine scaling approach to handle a wide class of constraints. We establish that CRNAS identifies a point satisfying ϵ-approximate second-order optimality conditions within O(ϵ−3/2) iterations. Finally, we compare CRNAS with MATLAB's optimization solver fmincon on three different test problems. These test problems all feature mixtures of heterogeneous populations, a problem setting that CRNAS is particularly well-suited for. Our numerical simulations show CRNAS has favorable performance, performing comparable if not better than fmincon in accuracy and computational cost for most of our examples.

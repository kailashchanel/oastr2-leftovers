# OASTR2 Leftovers

This is a codebase for a cosmology study primarily using Julia. This study aims to find the critical density (using percolation theory) that a galaxy cluster forms at. You can view the paper [here](insert link here (do not insert a link here until it is published)).

## Abstract
Numerical simulations and observations show that galaxies are not uniformly distributed. In cosmology, the largest known structures in the universe are galaxy filaments formed from the hierarchical clustering of galaxies due to gravitational forces. These consist of “walls” and “bridges” that connect larger superclusters with higher densities. Here, a graph theoretic approach is taken to model these structures as temporal networks. Using a method borrowed from statistical physics called percolation, cosmological graphs are reduced based on the valency of nodes to reveal the inner, most robust structural formation in time. By constraining the network, able to identify a threshold for physical features such as reachability, density, and connectivity, at which galactic matter allows for the formation of galaxy filaments in clusters can be identified. In this study, the infrastructure that may eventually be helpful in finding the threshold for physical features is created using the Julia programming language. 

## Codebase structure

This codebase is structured as a bunch of separate programs. There is no "one program" to run to run the entire analysis. Notes on what each program does can be found in a comment at the top each file.

Some code references a file called `GAMA_CZ5Unj.csv`. This file can be downloaded at the [GAMA (Galaxy And Mass Assembly survey) website](http://www.gama-survey.org/dr4/)

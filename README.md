# CLMassMod

## Code of Conduct

* We welcome all contributions from the community, and ask that
  behavior and interactions follow the LSST-DESC code of conduct:
  http://lsst-desc.org/sites/default/files/DESCProfessionalConductPolicies-2.pdf

## Contributing

* Alex Malz has put together a sketched out workflow for contributing to this repository:
https://docs.google.com/presentation/d/1UfrGJe2C-H9JVH6Gn2pklfbztF1tc7W4bVVMNszziiY/edit?usp=sharing 
* If you would like to be added as a collaborator, please contact one of the current admins (Camille Avestruz: avestruz@uchicago.edu, Alex Malz: aimalz@nyu.edu)

## Getting started: 
python setup.py develop

Currently, to use the Diemer+Kravtsov M-c relationship, CLMassMod also
requires a separate installation of Benedikt Diemer's code, colossus:
https://bitbucket.org/bdiemer/colossus

## Description
A galaxy cluster weak-lensing mass modeling and verification code

This repository documents the work done within the SPT-WL working
group and the LSST DESC Clusters working group to measure how weak
lensing mass estimates relate to the true 3D mass of a cluster. The
work focuses on emulating actual observations and mimicing procedures
used by weak lensing observers. Mock observations are based on
simulated lensing maps from n-body simulations.

## Set of to-do's for LSST-DESC sprint week's planning session

* Create a template for an ideal "config file" that the general user
 could modify

* Map out modules for the code (e.g. readers, observing strategies,
 M-c relations, (mis-)centering options and implementations, model
 profile to fit for mass (NFW, Diemer+, etc.), etc.

* Figure out how/where we can integrate with CCL and Tom M.'s
  cluster_toolkit package:
  https://github.com/tmcclintock/cluster_toolkit

* Map out steps to integrate SPHINX for documentation (Note - this
  should be done in parallel with refactoring):
  http://www.sphinx-doc.org/en/stable/contents.html

* Assign a list of mini-tasks so anyone can jump in

* Outline the process to jump in (submit an issue, reach out to
  "developers" involved for help, fork the repo, create a branch, make
  changes, merge into dev, etc.)




## Projects for the DESC Hack Day, 17 Feb 2017!

1.) Experiment with jointly fitting individual cluster shear profiles
(Mass, concentration, miscentering, 2-halo amplitude) jointly as a
population, ie a hierarchical model. Similar to [Lieu et al
2017](https://arxiv.org/abs/1701.00478), but we want to fit more
parameters per cluster and fit an order of magnitude more clusters
simultaneously. From a statistics & computation standpoint, can we
sample O(2k) parameters with an MCMC engine, perhaps by taking
advantage of the hierarchical structure of the problem?

2.) Daisuke Nagai & Erwin Lau, Nick Battaglia, and Stefan Hilbert &
the Magneticum team have offered lensing maps produced from hydro
simulations for analysis. Each simulation needs a reader that will
bring the lensing information to a common format for CLMASSMOD use. We
need to write that wrapper. We might also need to write an algorithm
that coverts 2D mass maps to shear maps.




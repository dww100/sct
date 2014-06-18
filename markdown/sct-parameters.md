# Format of SCT Parameter Files

SCT parameter files are encoded in [YAML](http://www.yaml.org/).
An example parameter file is provided [here](params.yml).
In this document we will describe all of the fields which can be included in the file and what they mean.
A level of familiarity with small angle scattering data and its analysis is assumed.
Many tutorials and review articles on such subjects are available, one such which connects the experimental data to the modelling analysis performed by SCT can be found [here](dx.doi.org/10.1016/j.ymeth.2011.01.004).

The first thing to note is that all distance units must be in Angstrom (Å).

## Curve Analysis

The first two sections in the example file are rg and rxs1. 
These sections contain parameters for the calculation of the for radius of gyration (Rg) and cross section (Rsx1) from scattering curves.

~~~~~~
---
rg:
    qmin: 0.0
    qmax: 0.2
    fitmin: 0.015
    fitmax: 0.028
rxs1:
    fitmin: 0.031
    fitmax: 0.051
~~~~~~

qmin and qmax are only used in plotting and represent the maximum and minimum values that will be plotted. They are ignored by all analysis scripts.

fitmin and fitmax are provided for both calculations, these are the minimum and maximum Q values used in the linear fits used to calculate the two quantities. These are used to analyse both the input experimental graphs and the theoretical curve generated from sphere modelling by SCT.

A third section can be added called rxs2 containing values for a second cross section at higher Q. 

## Sphere Modelling

The next section contains parameters used to generate sphere models (both dry for neutron and hydrated for x-ray comparisons)

~~~~~~
sphere:
    cutoff: 4
    boxside: 5.4454
~~~~~~

cutoff: Cut off number of atoms that must appear in a cube for a sphere to be included in the sphere model generated from an atomistic PDB model

boxside: Length of the side of the boxes used in the cubic grid used to create a sphere model from an atomistic PDB model

## Hydrating Sphere Models

Parameters used in the creation of the hydration layer around a sphere model (used for comparison to x-ray data)

~~~~~~
hydrate:
    positions: 27
    cutoff: 9
~~~~~~

positions: Number of positions on a box centred on each sphere from a dry model to place a preliminary hydration sphere (max = 27)

cutoff: Cut off number of spheres that must appear in a cube for a sphere to be included in the final hydration layer

## Theoretical Scattering Curve Calculation

Specify how the theoretical curve is generated from a sphere model.

~~~~~~
curve:
    qmax: 0.16
    npoints: 100
    radbins: 400
    smear: True
    wavelength: 6.0
    spread: 0.1
    divergence: 0.016
~~~~~~

qmax: Maximum Q value of the generated curve

npoints: Number of points to include in the curve

radbins: Number of bins used in creating separation histogram

smear: Will smearing be used in the creation of theoretical neutron curves (True or False)

wavelength: Smearing parameter

spread: Smearing parameter

divergence: Smearing parameter

## Curve Comparison

Range of Q values to use in the comparison of experimental and theoretical curves (i.e. to calculate the R factor or Chi^2 value).

~~~~~~
rfac:
    qmin: 0.015
    qmax: 0.16
~~~~~~

qmin: minimum Q value

qmax: maximum Q value


# Mountain-Building-and-Weathering-Efficiencies

## Overview
This repository contains contains the MATLAB code used for the climate modelling and weathering efficiency analysis presented within my research project. The code was originally developed by Rushby et al., (2018) and subsequently adapted for use in this research paper. 

## Files

### 'Updated_compute.m'
Main supporting function file which contains computational functions required by the main 'wrapper_bemrose.m' script. This file is not intended to be executed independently, instead, 'wrapper_bemrose.m' calls the relevant functions contained within it during execution.

### 'wrapper_bemrose.m'
Main script used to run analysis, including weathering efficiency, stellar and planetary parameters and heat flow.

### 'Weathering_Efficiency.m'
Main function for computing weathering efficiency values for any given radius, using an interpolated table from the values produced by W(R) in the main project.

## Requirements
MATLAB (Mathworks)

## How to Use
1) Open MATLAB
2) Add the repository folder to the MATLAB path
3) Run 'wrapper_bemrose.m'
4) Resulting outputs will correspond to analysis presented in this study

The code presented here is the version used for the final analysis in this research paper.

# PHD20 Firms and sectors in the macroeconomy
# Replication: Hopenhayn (1992) 
# February 2026

## 1. Packages
using Preferences, Parameters, Accessors, StaticArrays, Adapt, QuantEcon
using Base.Cartesian, LinearAlgebra, SparseArrays, LoopVectorization, Interpolations
using Distributions, Random, StatsBase, FastGaussQuadrature, Optim, Roots, Dierckx
using BenchmarkTools, AllocCheck, MAT

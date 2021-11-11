using LinearAlgebra
using DataFrames
using Random 
using Distributions 
using JLD 
using JuMP
using Ipopt 
#using SparseGrids
#using Latexify

Random.seed!(1)

T = 10; J = 10; K = 2; L = 2;
include("dataGeneratingFunc.jl");




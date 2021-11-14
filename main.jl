using Distributed # Introduce package for parallelization
addprocs(34) # Introduce CPU number
@everywhere cd(@__DIR__) 
@everywhere using LinearAlgebra
@everywhere using DataFrames
@everywhere using Random 
@everywhere using Distributions 
@everywhere using JLD 
@everywhere using JuMP
@everywhere using Ipopt 
#using SparseGrids
#using Latexify
@everywhere using StatsBase
@everywhere using Statistics
@everywhere Random.seed!(1)
@everywhere T = 10; @everywhere J = 10;@everywhere K = 2;@everywhere L = 2;


@sync @distributed for i = 81:97
    Random.seed!(1+i)
include("dataGenPar.jl");

save(join(["output/deneme",i,".jld"]), "charProduct_1",charProduct_1,"charProduct_2",charProduct_2,"df1",df1,"df2",df2);

end

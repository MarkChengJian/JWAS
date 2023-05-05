

using DataFrames              # package for working with data sets
using JWAS                    # package for Bayesian regression analyses, including BayesB and BayesCπ        
#using JWAS:misc               # utility functions
using JWAS.Datasets      
#using Plots                   # package for plotting 
using CSV
using DelimitedFiles

M = readdlm("genoJWAS.txt")
rowID = vec(readdlm("genoJWAS_ID.txt",String))

phenotypes = CSV.read("production_HIR_cycle5_7_11_12.dat", DataFrame,
    types=Dict(:ID => String), delim = ' ',header=true, missingstrings=["NA"] )
phenotypes= dropmissing(phenotypes, :QNurPenBatch)
phenotypes= dropmissing(phenotypes, :NurPenBatch)
phenotypes= dropmissing(phenotypes, :EntryAge)
phenotypes= dropmissing(phenotypes, :AMIR)

AMIR_model_equations = "AMIR = intercept + Batch + EntryAge + AMIR_0 + QNurPenBatch + SowID
                       nTrtsNur2_27 = intercept + Batch + EntryAge + NurPenBatch + SowID"
AMIR_R=[0.264518E-01 0
         0          1.1]
AMIR_model=build_model(AMIR_model_equations,AMIR_R);
set_covariate(AMIR_model,"EntryAge")
set_covariate(AMIR_model,"AMIR_0")

AMIR_G1=0.167989E-02
set_random(AMIR_model,"QNurPenBatch",AMIR_G1)
AMIR_G11=0.383042E-03
set_random(AMIR_model,"NurPenBatch",AMIR_G11)

AMIR_G2=[0.217102E-02  0
         0   0.193699E-01]
set_random(AMIR_model,"SowID",AMIR_G2)

AMIR_G3= [0.695238E-02   0
            0  0.147070]
@time add_genotypes(AMIR_model,M,AMIR_G3,header=false, rowID=rowID)

@time AMIR_outB=runMCMC(AMIR_model,phenotypes,methods="BayesB", missing_phenotypes=true,
    estimatePi=false,Pi=Dict([1.0; 1.0]=>0.001,[1.0; 0.0]=>0.001,[0.0; 1.0]=>0.001,[0.0; 0.0]=>0.998),
    estimateScale=true, chain_length=50000,burnin = 5000,output_heritability=true, 
    output_samples_frequency=100, outputEBV= true, output_folder = "AMIR_cNurTRT")
    
    out = GWAS("ChrInfo_map_1_7_Complete_JWAS_QC.txt", 
    "AMIR_cNurTRT/MCMC_samples_marker_effects_geno_AMIR.txt",
    "AMIR_cNurTRT/MCMC_samples_marker_effects_geno_nTrtsNur2_27.txt",
    AMIR_model,header=true,window_size="1 Mb",
    GWAS=false,sliding_window=false,genetic_correlation=true)
    
CSV.write("MCMC_BayesB_1Mb_window_AMIR_cNurTRT.txt", out[1])



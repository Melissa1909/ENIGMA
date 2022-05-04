clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));

% import cortical data
CT = readmatrix('T_lhrh_thickness_term-preterm.csv');
GYR = readmatrix('T_lhrh_gyrification_term-preterm.csv');
genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");


outpath = ('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\ENIGMA\matlab\results\cellTypes');


%% fetch celltype-specific genes based on Li 2021
cellGenes = GiveMeCellGenes();

%% perform correlation with phenotype
cellGenes=fopen(fullfile(outpath, "Astro.txt"));

corr_coeffs_CT = NaN(numGenes,6);

for g = 1:numGenes
    fprintf('Performing correlation for gene %d\n', g);
    geneName = astro(g);
    y=strcmp(geneName,genes.Var1);
    geneHere = genes(y,2:69);
    geneHere = table2array(geneHere);
    [rho,pval] = corr(CT',geneHere','Type','Spearman');
    pspin = spin_test(CT',geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
    
    %save output
    corr_coeffs_CT(g,2) = rho;
    corr_coeffs_CT(g,3) = pval;
    corr_coeffs_CT(g,5) = pspin;
end


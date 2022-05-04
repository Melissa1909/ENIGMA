clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));


genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

outpath = ('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\ENIGMA\matlab\results\cellTypes');

measure = "CT";
celltype = "Astro";
celltypes=["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"];


%% fetch celltype-specific genes based on Li 2021
cellGenes = GiveMeCellGenes(celltype);

%% perform correlation with phenotype
corr_coeffs = GiveMeCorrelation(measure,cellGenes,genes_all);

%% multiple comparison correction
corr_coeffs = array2table(corr_coeffs);
corr_coeffs = renamevars(corr_coeffs,["corr_coeffs1","corr_coeffs2","corr_coeffs3","corr_coeffs4","corr_coeffs5","corr_coeffs6"], ...
["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
corr_coeffs.gene_symbol = cellGenes;

% adjust p for FDR
[p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p, 0.05, false);
corr_coeffs.p_fdr = p_fdr';

% adjust p_spin for FDR
[pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p_spin, 0.05, false);
corr_coeffs.p_spin_fdr = pspin_fdr';

%% save 
fileNameOut = sprintf('correlation_%s_%s.csv',measure,celltype);
fileNameOut = fullfile(outpath,fileNameOut);
writetable(corr_coeffs,fileNameOut);


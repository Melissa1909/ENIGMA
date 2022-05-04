clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));


% Import data
genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));
CT = readmatrix('T_lhrh_thickness_term-preterm.csv');
GYR = readmatrix('T_lhrh_gyrification_term-preterm.csv');

%% correlation
numGenes = 10; %lower the number for debugging and test purposes
%numGenes = height(genes);

corr_coeffs_CT = NaN(numGenes,6);


for g = 1:numGenes
    fprintf('Performing correlation for gene %d\n', g);
    geneHere = genes(g,:);
    geneHere = table2array(geneHere);
    [rho,pval] = corr(CT',geneHere','Type','Spearman');
    pspin = spin_test(CT',geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
    
    %save output
    corr_coeffs_CT(g,2) = rho;
    corr_coeffs_CT(g,3) = pval;
    corr_coeffs_CT(g,5) = pspin;
end


%% MCC old

% adjust p for FDR
%corr_coeffs_CT.p_fdr = mafdr(corr_coeffs_CT.p);

% adjust p_spin for FDR
%corr_coeffs_CT.p_spin_fdr = mafdr(corr_coeffs_CT.p_spin);

%% MCC
corr_coeffs_CT = array2table(corr_coeffs_CT);
corr_coeffs_CT = renamevars(corr_coeffs_CT,["corr_coeffs_CT1","corr_coeffs_CT2","corr_coeffs_CT3","corr_coeffs_CT4","corr_coeffs_CT5","corr_coeffs_CT6"], ...
["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
corr_coeffs_CT.gene_symbol = genelabels(1:numGenes);

% adjust p for FDR
[p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs_CT.p, 0.05, false);
corr_coeffs_CT.p_fdr = p_fdr';

% adjust p_spin for FDR
[pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs_CT.p_spin, 0.05, false);
corr_coeffs_CT.p_spin_fdr = pspin_fdr';

% order data
corr_coeffs_CT_sorted = sortrows(corr_coeffs_CT,"p_spin_fdr");


clear; clc;
startup;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\sw\BrainSpace-0.1.2'))

% Import data
genes_all = readtable('expression_brainorder_lh-mirror.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

CT_all = readtable('thickness_tvals_multipleRegression_combatCorr.csv');
CT = CT_all.tval;
GYR_all = readtable('gyrification_tvals_multipleRegression_combatCorr.csv');
GYR = GYR_all.tval;

%% correlation
numGenes = 10; %lower the number for debugging and test purposes
%numGenes = height(genes);

corr_coeffs_CT = NaN(numGenes,6);

measures = ["CT", "GYR"];
for measureNum=1:length(measures)
    measure=measures(measureNum);

    for g = 1:numGenes
        fprintf('Performing correlation for gene %d\n', g);
        geneHere = genes(g,:);
        geneHere = table2array(geneHere);
        [rho,pval] = corr(CT,geneHere','Type','Spearman');
        pspin = spin_test(CT,geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
        
        %save output
        corr_coeffs_CT(g,2) = rho;
        corr_coeffs_CT(g,3) = pval;
        corr_coeffs_CT(g,5) = pspin;
    end
   
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


% 
% AUTHOR: MELISSA THALHAMMER
%
% This script uses Gyr and CTh group difference between preterm and
% term-born subjects and correlates these Tvals with gene expression.
%
% TO DO: clean the script so that paths go to analysis2/data


clear; clc;
startup;
working_dir='C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/';
addpath(genpath(working_dir));
addpath(genpath(fullfile(working_dir, 'data_melissa')));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\sw\BrainSpace-0.1.2'));

% Import data
genes_all = readtable('expression_brainorder_lh-mirror.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

CT_all = readtable('thickness_tvals_multipleRegression_combatCorr.csv');
CT = CT_all.tval;
GYR_all = readtable('gyrification_tvals_multipleRegression_combatCorr.csv');
GYR = GYR_all.tval;

%% correlation
numGenes = 100; %lower the number for debugging and test purposes
%numGenes = height(genes);

measures = [CT, GYR];
for measureNum=1:width(measures)
    measure=measures(:,measureNum);
    corr_coeffs = NaN(numGenes,6);

    for g = 1:numGenes
        fprintf('Performing correlation for gene %d\n', g);
        geneHere = genes(g,:);
        geneHere = table2array(geneHere);
        [rho,pval] = corr(measure,geneHere','Type','Spearman');
        [pspin, r_dist] = spin_test(measure,geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');

        %save output
        corr_coeffs(g,2) = rho;
        corr_coeffs(g,3) = pval;
        corr_coeffs(g,5) = pspin;
    end
   
    corr_coeffs = array2table(corr_coeffs);
    corr_coeffs = renamevars(corr_coeffs,["corr_coeffs1","corr_coeffs2","corr_coeffs3","corr_coeffs4","corr_coeffs5","corr_coeffs6"], ...
    ["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
    corr_coeffs.gene_symbol = genelabels(1:numGenes);
    
    % adjust p for FDR
    [p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p, 0.05, false);
    corr_coeffs.p_fdr = p_fdr';
    
    % adjust p_spin for FDR
    [pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p_spin, 0.05, false);
    corr_coeffs.p_spin_fdr = pspin_fdr';
    
    % order data
    corr_coeffs_sorted = sortrows(corr_coeffs,"p_spin_fdr");

     if measureNum == 1
        m = "CT";
    elseif measureNum == 2
        m = "GYR";
    end
    %writetable(corr_coeffs_sorted, fullfile(working_dir, 'results\p_spin_corr\',m));


end

%% DONE!
fprintf("Done!");
function corr_coeffs = singleWorkflow(measure,celltype,genes_all,outpath)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% Performs:
%       1. Gene fetching for specified cell type,
%       2. ADDED: average gene expression of one cell type as zscore
%       3. Spearman correlation of cortical measure and gene expression for specified cell type,
%       4. FDR-correction of pval and pspin
%
%---INPUTS:
% * measure: "CT" or "GYR".
% * celltype: one of the following: ["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"].
% * genes_all: AHBA gene expression data, preprocessed with abagen. 
%               Format should be (numGenes, DK-Rois).
% * outpath: path where data should be saved.
%
% 
%---OUTPUTS:
% * corr_coeffs: contains Spearman rho, parametric pval, pspin. 

%----------------------------------------------------------------------


%% fetch celltype-specific genes based on Li 2021
cellGenes = GiveMeCellGenes(celltype);
% for test purposes:
%cellGenes = ["A1BG", "ACAA2", "thisgenedoesnotexist"]';

%% perform correlation with phenotype
corr_coeffs = GiveMeCorrelation(measure,cellGenes,genes_all);

%% multiple comparison correction
corr_coeffs = array2table(corr_coeffs);
corr_coeffs = renamevars(corr_coeffs,["corr_coeffs1","corr_coeffs2","corr_coeffs3","corr_coeffs4","corr_coeffs5","corr_coeffs6"], ...
["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
corr_coeffs.gene_symbol = cellGenes';

% adjust p for FDR
[p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p, 0.05, false);
corr_coeffs.p_fdr = p_fdr';

% adjust p_spin for FDR
[pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p_spin, 0.05, false);
corr_coeffs.p_spin_fdr = pspin_fdr';

% save 
fileNameOut = sprintf('correlation_%s_%s.csv',measure,celltype);
fileNameOut = fullfile(outpath,fileNameOut);
corr_coeffs_sorted = sortrows(corr_coeffs,"p_spin_fdr");
writetable(corr_coeffs_sorted,fileNameOut);

%% dropna
corr_coeffs(any(ismissing(corr_coeffs), 2),:) = [];
fprintf('For celltype %s, only %i were found in the AHBA...\n', celltype, size(corr_coeffs,1))
% adjust p for FDR
[p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p, 0.05, false);
corr_coeffs.p_fdr = p_fdr';

% adjust p_spin for FDR
[pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs.p_spin, 0.05, false);
corr_coeffs.p_spin_fdr = pspin_fdr';

% save droppedna version
fileNameOut = sprintf('correlation_%s_%s_dropna.csv',measure,celltype);
fileNameOut = fullfile(outpath,fileNameOut);
corr_coeffs_sorted = sortrows(corr_coeffs,"p_spin_fdr");
writetable(corr_coeffs_sorted,fileNameOut);



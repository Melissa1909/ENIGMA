function corr_coeffs = singleWorkflow_DiBiase(measure,celltype,genes_all,outDir)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% Performs:
%       1. Gene fetching for specified cell type,
%       2. ADDED: average gene expression of one cell type as zscore
%       3. Spearman correlation of cortical measure and gene expression for specified cell type,
%       
%
%---INPUTS:
% * measure: "CT" or "GYR".
% * celltype: one of the following: ["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"].
% * genes_all: AHBA gene expression data, preprocessed with abagen. 
%               Format should be (numGenes, DK-Rois).
% * outDir: path where data should be saved.
%
% 
%---OUTPUTS:
% * corr_coeffs: contains Spearman rho, parametric pval, pspin. 

%----------------------------------------------------------------------


%% fetch celltype-specific genes based on Li, Seidlitz 2021
if celltype=="Subplate"
    cellGenes = readtable("C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\analysis\results\hypothesis-driven\pathophysiology\subplate_neurons\tested_genes.txt", 'Format','%s','ReadVariableNames', false);
    cellGenes = rows2vars(cellGenes);
    cellGenes = table2array(cellGenes(:,[2:width(cellGenes)]));
    % remove empty elements
    vb = cellfun(@isempty, cellGenes);
    cellGenes = cellGenes(~vb);
    % remove duplicates
    cellGenes = unique(cellGenes,'stable');
    cellGenes = convertCharsToStrings(cellGenes);
else
    cellGenes = GiveMeCellGenes(celltype, outDir);
end

% for test purposes:
%cellGenes = ["A1BG", "ACAA2", "thisgenedoesnotexist"]';

%% Mean expression of each cell type-specific gene set, z-scored
expression = GiveMeCellExpression(cellGenes, genes_all);

%% perform correlation with phenotype
corr_coeffs = GiveMeCellCorrelation(measure,expression);

corr_coeffs = array2table(corr_coeffs);
corr_coeffs = renamevars(corr_coeffs,["corr_coeffs1","corr_coeffs2","corr_coeffs3"], ...
["Spearman_rho","p","p_spin"]);


% save 
fileNameOut = sprintf('correlation_cellTypeMean_%s_%s.csv',measure,celltype);
fileNameOut = fullfile(outDir,fileNameOut);
writetable(corr_coeffs,fileNameOut);





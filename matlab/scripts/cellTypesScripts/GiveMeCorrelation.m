function corr_coeffs = GiveMeCorrelation(measure, cellGenes, genes_all)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% Performs Spearman correlation: parametric one and pspin based on ENIGMA
% toolbox. 
%
%---INPUTS:
% * measure: "CT" or "GYR".
% * cellGenes: to be loaded using `GiveMeCellGenes.m` for one celltype.
% * genes_all: AHBA gene expression data, preprocessed with abagen. 
%               Format should be (numGenes, DK-Rois).
% 
%---OUTPUTS:
% * corr_coeffs: contains Spearman rho, parametric pval, pspin. 

%----------------------------------------------------------------------
% Define measure
if measure == "CT"
    corticalData = readmatrix('T_lhrh_thickness_term-preterm.csv');
elseif measure == "GYR" 
    corticalData = readmatrix('T_lhrh_gyrification_term-preterm.csv');
else
    fprintf('Please specify a cortical measure...')
end

if ~(isequal(size(corticalData,2),68))
    warning('ERROR WHEN READING DATA! RESHAPING ARRAY...');
    corticalData = corticalData(2:69);
end


fprintf('Starting calculations for measure %s\n', measure);

% output array
numGenes = length(cellGenes);
corr_coeffs = NaN(numGenes,6);


%----------------------------------------------------------------
for g = 1:numGenes
    fprintf('Performing correlation for gene %d\n', g);
    geneName = cellGenes(g);
    y=strcmp(geneName,genes_all.Var1);
    geneHere = genes_all(y,2:69);
    geneHere = table2array(geneHere);
    if any(geneHere)
        [rho,pval] = corr(corticalData',geneHere','Type','Spearman');
        pspin = spin_test(corticalData',geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
        
        %save output
        corr_coeffs(g,2) = rho;
        corr_coeffs(g,3) = pval;
        corr_coeffs(g,5) = pspin;
    else
        corr_coeffs(g,:) = nan;
    end
    
end
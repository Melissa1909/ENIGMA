function expressionAllMeanZ = GiveMeCellExpression(cellGenes, genes_all)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% Based on the idea presented in Di Biase 2022: mean expression of each
% cell type-specific gene set is determined within each DK-ROI and
% normalized (z-scored) by mean-expression across the entire brain. 
%
%---INPUTS:
% * cellGenes: to be loaded using `GiveMeCellGenes.m` for one celltype.
% * genes_all: AHBA gene expression data, preprocessed with abagen. 
%               Format should be (numGenes, DK-Rois).
% 
%---OUTPUTS:
% * expression: cell-type specific expression values of all genes relevant
%               for cell type, summarized by mean. 

%----------------------------------------------------------------------

numGenes = length(cellGenes);

%% calculate mean expression of each cell-type gene set within each DK-ROI
expressionAll = NaN(numGenes,68);

for g = 1:numGenes
    % get expression value of each relevant gene for the cell type
    geneName = cellGenes(g);
    y=strcmp(geneName,genes_all.Var1);

    % check whether gene name is in AHBA genes_all
    if any(y)
        geneHere = genes_all(y,2:69);
        geneHere = table2array(geneHere);
    else
        geneHere = nan;
    end

    % store expression values for all genes of relevant cell type
    expressionAll(g,:) = geneHere;
end

% dropna rows
expressionAll(any(ismissing(expressionAll), 2),:) = [];

% get mean expression for gene set for each DK-ROI
expressionAllMean = mean(expressionAll, 1);

%% normalization using z-score --> IS THAT CORRECT?
%warning('This is my interpretation of the Di Biase paper...');
%meanSample = mean(expressionAllMean,2); %mean expression across the brain
%stdSample = std(expressionAllMean,0,2);

%expressionAllMeanZ = NaN(1,68);
%for row=1:68
%    expressionAllMeanZ(1,row) = (expressionAllMean(1,row) - meanSample) / stdSample;
%end

expressionAllMeanZ = expressionAllMean;





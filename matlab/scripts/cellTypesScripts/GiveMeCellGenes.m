function cellGenes = GiveMeCellGenes(c, outDir)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% determines all genes relevant for one cell class
% data and idea from Li, Seidlitz 2021
%
%---INPUTS:
% * c: cell type of list ["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro",
%       "Oligo", "OPC", "Per"].
% * outDir: path where data should be saved.

%----------------------------------------------------------------------
% Read in data and check inputs
cell_raw = readtable('cell_types_Li-Seidlitz2021_melissa.csv', 'ReadVariableNames',0);

if nargin < 1
    error('You must specify a cell type from the given list...');
end

%---------------------------------------------------------------------------
cellRow = [];
    
% mark all lines that belong to current celltype
for rows=1:height(cell_raw)
    currentLine=cell_raw(rows,:);
    if (strcmp(currentLine.Var4,c))  
        cellRow(rows,1)=1;
        %cellRow = logical(cellRow);
    end
end
cellRow = logical(cellRow);

% store all genes that belong to the current celltype
helperGenes={};
helperGenes=cell_raw(cellRow,[5:713]);
helperGenes=table2array(helperGenes);

cellGenes=strings(1,1);
geneIndex=1;
for a=1:size(helperGenes,1)
    for b=1:size(helperGenes,2)
        cellGenes(1,geneIndex)=helperGenes(a,b);
        geneIndex=geneIndex+1;
    end
end

% remove empty elements
vb = cellfun(@isempty, cellGenes);
cellGenes = cellGenes(~vb);
% remove duplicates
cellGenes = unique(cellGenes,'stable');
     
% save output
outname=fullfile(outDir,c);
writematrix(cellGenes,outname);
     

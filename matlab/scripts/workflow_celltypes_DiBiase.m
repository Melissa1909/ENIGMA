% This version works in analogy to the paper Di Biase 2022. 
% 1. A cell-type-specific gene set is determined based on some single-cell
%    seq data. Here: Li, Seidlitz 2021 have provided a link to a
%    publication-based summary of relevant genes for each of the 7 cell
%    types (Di Biase seem to have extended that one, also downloadable,
%    which would be one possible update of this script).
%       -> GiveMeCellGenes.m
% 
% 2. Mean expression of each cell type-specific gene set, z-scored, is
%    determined. 
%       -> GiveMeCellExpression.m
%
% 3. Correlation for both measures with each of the 7 cell-type-specific
%    map is performed, respectively.
%


clear; clc;
projectDir = ('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\Melissa1909_ENIGMA\matlab');
addpath(genpath(projectDir));
outDir = fullfile(projectDir, 'results/cellTypes/');

startup;

genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

measures = ["CT", "GYR"];
%measures = ["CT"];
celltype = "Astro";
celltypes=["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"];

% test
%singleWorkflow_DiBiase(measures,celltype,genes_all,outDir);

for measureNum=1:length(measures)
    measure=measures(measureNum);

    for cellNum=1:length(celltypes)
        celltype = celltypes(cellNum);
        fprintf('STARTING WITH NEW CELLTYPE: %s\n', celltype);

        singleWorkflow_DiBiase(measure,celltype,genes_all,outDir);
    end
end

clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
%addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));
startup;

genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

outpath = ('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\ENIGMA\matlab\results\cellTypes');

%measures = ["CT", "GYR"];
measures = ["CT"];
celltype = "Astro";
celltypes=["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"];

% test
%singleWorkflow(measures,celltype,genes_all,outpath);

for measureNum=1:length(measures)
    measure=measures(measureNum);

    for cellNum=1:length(celltypes)
        celltype = celltypes(cellNum);
        fprintf('STARTING WITH NEW CELLTYPE: %s\n', celltype);

        singleWorkflow(measure,celltype,genes_all,outpath);
    end
end

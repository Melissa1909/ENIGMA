clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/data_melissa'));

% import cortical data
CT = readmatrix('T_lhrh_thickness_term-preterm.csv');
GYR = readmatrix('T_lhrh_gyrification_term-preterm.csv');
genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");

cell_raw = readtable('cell_types_Li-Seidlitz2021_melissa.csv', 'ReadVariableNames',0);
% import cell type genes
%astro=readtable('astro_test.csv');
astro=["SLC1A2","RNF219-AS1","SLC1A3","COL5A3","ATP1A2"];
%% cell_raw = readtable('cell_types_Li-Seidlitz2021_melissa.csv');
astroRow=[];
endoRow=[];
neuroexRow=[];
neuroinRow=[];
microRow=[];
oligoRow=[];
opcRow=[];
perRow=[];


for rows = 1:height(cell_raw)
    currentLine=cell_raw(rows,:);
     if (strcmp(currentLine.Class,'Astro'))
        astroRow(rows,1)=1;
        astroRow = logical(astroRow);
     elseif (strcmp(currentLine.Class,'Endo'))
         endoRow(rows,1)=1;
         endoRow = logical(endoRow);
     elseif (strcmp(currentLine.Class,'Neuro-Ex'))
         neuroexRow(rows,1)=1;
         neuroexRow = logical(neuroexRow);
     elseif (strcmp(currentLine.Class,'Neuro-In'))
         neuroinRow(rows,1)=1;
         neuroinRow = logical(neuroinRow);
     elseif (strcmp(currentLine.Class,'Micro'))
         MicroRow(rows,1)=1;
         MicroRow = logical(MicroRow);
     elseif (strcmp(currentLine.Class,'Oligo'))
         oligoRow(rows,1)=1;
         oligoRow = logical(oligoRow);
     elseif (strcmp(currentLine.Class,'OPC'))
         opcRow(rows,1)=1;
         opcRow = logical(opcRow);
     elseif (strcmp(currentLine.Class,'Per'))
         perRow(rows,1)=1;
         perRow = logical(perRow);
    end
end
astro=cell_raw(astroRow,[5:713]);
astroGenes=[];
for g=1:length(astro)
    for h=1:height(astro)
        astroGenes(rows)=g
    end
end



endo=cell_raw(endoRow,:);
ex=cell_raw(neuroexRow,:);
in=cell_raw(neuroinRow,:);
micro=cell_raw(microRow,:);
oligo=cell_raw(oligoRow,:);
opc=cell_raw(opcRow,:);
per=cell_raw(perRow,:);
%% loop
celltypes=["Astro", "Endo"];
cellGenes=[];

for c=celltypes
    %cellRow=[];
    for rows = 1:height(cell_raw)
        currentLine=cell_raw(rows,:);
        if (strcmp(currentLine.Var4,c))  
            cellRow(rows,1)=1;
            cellRow = logical(cellRow);
          %cellGenes=cell_raw(rows,[5:713]);
          %cellGenes = rmmissing(cellGenes,2);

        end
        cellGenes=cell_raw(cellRow,[5:713]);
        cellGenes = rmmissing(cellGenes,2);
    end
    %vector=[];
    %for i=1:size(cellGenes,1)
       % vector=[vector, cellGenes(i,:)];
    %end
    cellGenes2=reshape(cellGenes, [1,288])
end



%% correlation


%numGenes = 10; %lower the number for debugging and test purposes
numGenes = height(astro);

corr_coeffs_CT = NaN(numGenes,6);

for g = 1:numGenes
    fprintf('Performing correlation for gene %d\n', g);
    geneName = astro(g);
    y=strcmp(geneName,genes.Var1);
    geneHere = genes(y,2:69);
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



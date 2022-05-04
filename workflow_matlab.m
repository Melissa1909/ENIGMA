%addpath(genpath('path-to-enigma/ENIGMA/matlab/'));
clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
    
%% Gene expression data
% Fetch gene expression data
genes_all = fetch_ahba();
genes_all = rows2vars(genes_all);        
genes = genes_all(2:12669,2:69);

% Obtain region labels
%reglabels = genes.label;

% Obtain gene labels
%genelabels = genes.Properties.VariableNames(2:end);

%% Simulate Cortical Thickness Data for 68 regions and 200 subs
raw_ct = (rand(68,1))+2;
%raw_ct = raw_ct';
test_data = (rand(68,1))+3;

%% Correlation
corr_coeffs_CT = NaN(2,6);
test = corr(geneHere',raw_ct);

rows = height(genes);
for g = 2:3
    fprintf('Performing correlation for gene %d\n', g);
    geneHere = genes(g,:);
    geneHere = table2array(geneHere);
    geneHere = cell2mat(geneHere);
    X = raw_ct';
    Y = geneHere';
    [rho,pval] = corr(X,Y,'Type','Spearman');
    pspin = spin_test(raw_ct',geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
    
    %save output
    corr_coeffs_CT(g,2) = rho;
    corr_coeffs_CT(g,3) = pval;
    corr_coeffs_CT(g,5) = pspin;
end

%% MCC
corr_coeffs_CT = array2table(corr_coeffs_CT);
corr_coeffs_CT = renamevars(corr_coeffs_CT,["corr_coeffs_CT1","corr_coeffs_CT2","corr_coeffs_CT3","corr_coeffs_CT4","corr_coeffs_CT5","corr_coeffs_CT6"], ...
    ["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
%corr_coeffs_CT.gene_symbol = genelabels;

% adjust p for FDR
corr_coeffs_CT.p_fdr = mafdr(corr_coeffs_CT.p);

% adjust p_spin for FDR
corr_coeffs_CT.p_spin_fdr = mafdr(corr_coeffs_CT.p_spin);



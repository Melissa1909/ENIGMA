%addpath(genpath('path-to-enigma/ENIGMA/matlab/'));
clear; clc;
addpath(genpath('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression/ENIGMA/matlab/'));
    
%% Gene expression data
genes_all = readtable('expression_brainorder.csv','ReadVariableNames',0);
genelabels = genes_all.Var1;
genes = removevars(genes_all, "Var1");


%% Simulate Cortical Thickness Data for 68 regions
%raw_ct = (rand(68,1))+2;
raw_ct = [2.58805912205776;2.47415092069538;2.91671074814889;2.70854757827065;2.73262913297386;2.01698900301724;2.11276871049553;2.66646281648323;2.41444090581014;2.02080723998185;2.04091515917065;2.01846176337805;2.54068235120194;2.98425707100468;2.07604808494799;2.88011775234679;2.83816330333713;2.32463650216453;2.74662763295216;2.25655974079322;2.90183144477258;2.22360150324691;2.98971228734974;2.05102334396157;2.76607876404913;2.70086328634064;2.20483752327492;2.49138462181976;2.64045953540602;2.10324598242545;2.63940741103014;2.50666531421732;2.64030120700556;2.20503425655978;2.57767854951135;2.21498410838450;2.80366321980188;2.96376106621246;2.04047398786189;2.70892046478552;2.98341210942531;2.37947296931904;2.14420901734182;2.62838728785446;2.83654282640552;2.96933843689782;2.22493682200663;2.49265154232225;2.57427236429738;2.43024895189780;2.57240342057358;2.96272977288615;2.25077821778144;2.70811834343147;2.21925017057708;2.02388433108541;2.23236613457015;2.83288178787764;2.87043287640347;2.33179605075937;2.21401196526758;2.36595076912635;2.41593958701836;2.89946720368525;2.93131599653450;2.02013578431943;2.38285745191099;2.27187183970819];

%% Correlation
numGenes = 20;

corr_coeffs_CT = NaN(numGenes,6);

rows = height(genes);
for g = 1:numGenes
    fprintf('Performing correlation for gene %d\n', g);
    geneHere = genes(g,:);
    geneHere = table2array(geneHere);
    [rho,pval] = corr(raw_ct,geneHere','Type','Spearman');
    pspin = spin_test(raw_ct,geneHere','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
    
    %save output
    corr_coeffs_CT(g,2) = rho;
    corr_coeffs_CT(g,3) = pval;
    corr_coeffs_CT(g,5) = pspin;
end

%% MCC
corr_coeffs_CT = array2table(corr_coeffs_CT);
corr_coeffs_CT = renamevars(corr_coeffs_CT,["corr_coeffs_CT1","corr_coeffs_CT2","corr_coeffs_CT3","corr_coeffs_CT4","corr_coeffs_CT5","corr_coeffs_CT6"], ...
    ["gene_symbol","Spearman_rho","p","p_fdr","p_spin","p_spin_fdr"]);
corr_coeffs_CT.gene_symbol = genelabels(1:numGenes);

% adjust p for FDR
%corr_coeffs_CT.p_fdr = mafdr(corr_coeffs_CT.p);
[p_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs_CT.p, 0.05, false);
corr_coeffs_CT.p_fdr = p_fdr';

% adjust p_spin for FDR
%corr_coeffs_CT.p_spin_fdr = mafdr(corr_coeffs_CT.p_spin);
[pspin_fdr, c_alpha, h, extra] = fdr_BH(corr_coeffs_CT.p_spin, 0.05, false);
corr_coeffs_CT.p_spin_fdr = pspin_fdr';





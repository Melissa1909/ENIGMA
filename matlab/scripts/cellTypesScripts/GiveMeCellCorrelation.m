function corr_coeffs = GiveMeCellCorrelation(measure,expression)
%----------------------------------------------------------------------
% AUTHOR: MELISSA THALHAMMER
%
% Performs Spearman correlation: parametric one and pspin based on ENIGMA
% toolbox. 
%
%---INPUTS:
% * measure: "CT" or "GYR".
% * expression: expression of cell-type specific gene set. Obtained from
%               GiveMeCellExpression.m
% 
%---OUTPUTS:
% * corr_coeffs: contains Spearman rho, parametric pval, pspin. 

%----------------------------------------------------------------------

% Define measure
if measure == "CT"
    %corticalData = readmatrix('T_lhrh_thickness_term-preterm.csv');
    corticalData = readtable('thickness_tvals_multipleRegression_combatCorr');
    corticalData = corticalData.tval';
elseif measure == "GYR" 
    corticalData = readtable('gyrification_tvals_multipleRegression_combatCorr.csv');
    corticalData = corticalData.tval';
else
    fprintf('Please specify a cortical measure...')
end

if ~(isequal(size(corticalData,2),68))
    warning('ERROR WHEN READING DATA! RESHAPING ARRAY...');
    corticalData = corticalData(2:69);
end


fprintf('Starting calculations for measure %s\n', measure);

% output array
corr_coeffs = NaN(1,3);


%----------------------------------------------------------------

[rho,pval] = corr(corticalData',expression','Type','Spearman');
pspin = spin_test(corticalData',expression','surface_name', 'fsa5', 'parcellation_name', 'aparc', 'n_rot', 1000,'type', 'spearman');
        
%save output
corr_coeffs(1,1) = rho;
corr_coeffs(1,2) = pval;
corr_coeffs(1,3) = pspin;
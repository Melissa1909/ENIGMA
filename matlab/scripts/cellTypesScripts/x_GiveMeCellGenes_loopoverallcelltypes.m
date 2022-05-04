function cellGenes = GiveMeCellGenes()
%----------------------------------------------------------------------
% determines all genes relevant for one cell class
% data and idea from Li, Seidlitz 2021


cell_raw = readtable('cell_types_Li-Seidlitz2021_melissa.csv', 'ReadVariableNames',0);
celltypes=["Astro", "Endo", "Neuro-Ex", "Neuro-In", "Micro", "Oligo", "OPC", "Per"];
outpath = ('C:\Users\Acer\Documents\Studium\PhD\01_MA_preterm_gene-expression\ENIGMA\matlab\results\cellTypes');

%---------------------------------------------------------------------------

for c=celltypes
    cellRow = [];
    
    % mark all lines that belong to current celltype
    for rows=1:height(cell_raw)
        currentLine=cell_raw(rows,:);
        if (strcmp(currentLine.Var4,c))  
            cellRow(rows,1)=1;
            cellRow = logical(cellRow);
        end
    end

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
     outname=fullfile(outpath,c);
     writematrix(cellGenes,outname);
     
end
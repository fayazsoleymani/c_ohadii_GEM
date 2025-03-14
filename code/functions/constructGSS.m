function GSS= constructGSS(GSSfilePath)
    % description: This function construct a GSS struction in order to be
    % used for compartmentalization
    % inputs
    %   GSSfilePath     path to the GSS file, in which the probabilites are
    %                   divided by the most probable compartment, so that
    %                   the most probable location of the gene is equal to
    %                   one. 
    % outputs
    %   GSS             the GSS struct.

    GSSdata= readtable(GSSfilePath);
    GSS= struct();
    if ~isfield(GSS, 'compartmets')
        GSS.compartments= GSSdata.Properties.VariableNames(2:end);
    end
    if ~isfield(GSS, 'scores')
        GSS.scores= table2array(GSSdata(:, 2:end));
    end
    if ~isfield(GSS, 'genes')
        GSS.genes= table2cell(GSSdata(:, 1));
    end
    if ~isfield(GSS, 'geneComp')
        GSS.geneComp= zeros(size(GSS.genes, 1), 1);
    end
    for geneIdx= 1:numel(GSS.genes)
        comp= find(GSS.scores(geneIdx, :)==1);
        if length(comp) > 1
            comp=1;
        end
        GSS.geneComp(geneIdx) = comp;
    end
end
function model= addBiomassRxn(model, biomassFile)
    % function to add the biomass reaction based on the input file
    % input
    %       model           input model
    %       biomassFile     the file with all biomass coefficients, the
    %                       column names are the biomass reaction names
    % output
    %       model           output model
    
    biomassTable= readtable(biomassFile, 'Delimiter', '\t', 'FileType', ...
        'text','VariableNamingRule', 'preserve');
    biomassNames= biomassTable.Properties.VariableNames(2:end);
    for biomassIndex= 1: length(biomassNames)
        biomassName= biomassNames{biomassIndex};
        coefsDict= dictionary();

        for index= 1:numel(biomassTable.metName)
            metName= biomassTable.metName{index};
            metIndex= find(strcmp(model.metNames, metName));
            coefsDict(model.mets(metIndex))= -biomassTable.(biomassName)(index);
        end

        model = addReaction(model, biomassName, ...
            'metaboliteList', keys(coefsDict), ...
            'stoichCoeffList', values(coefsDict),...
            'reversible', false, 'lowerBound', 0, 'upperBound', 1000);
    end
end
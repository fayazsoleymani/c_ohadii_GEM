function model= addRxnsAnnotations(model, annotationFile)
    % description: add the reaction annotations to the model
    % inputs
    %       model           input model
    %       annotationFile  file containing the annotations, e.g. EC-number
    %       , KEGG reaction, etc.
    % output
    %       model           output model

    mappingTable= readtable(annotationFile);

    for i= 1:size(mappingTable, 1)
        rxnIndex= find(strcmp(mappingTable.ObjectID{i}, model.rxns));
        if ~isempty(mappingTable.EC_Number{i})
            try
                model = addMIRIAMAnnotations(model, model.rxns(rxnIndex), ...
                    'ec-code' , mappingTable.EC_Number(i));
            catch
                continue
            end
        end
        if ~isempty(mappingTable.KEGG_reaction{i})
            model= addMIRIAMAnnotations(model, model.rxns(rxnIndex),...
                'kegg.reaction', mappingTable.KEGG_reaction(i));
        end
    end

end
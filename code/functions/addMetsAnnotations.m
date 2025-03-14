function model= addMetsAnnotations(model, annotationFile)
    % description: add the metabolite annotations to the model
    % inputs
    %       model           input model
    %       annotationFile  file containing the annotations, e.g. ChEBI,
    %                       BIGG, ChemSpider, MetaNetX, PubChem, SEED,
    %                       InChI, InChIKey, KEGG, 
    % output
    %       model           output model


    clearedCompsModel= clearComps(model, false, true, false);
    clearedMets= clearedCompsModel.mets;
    mappingTable= readtable(annotationFile);

    for i= 1:size(mappingTable, 1)
        metIndices= find(strcmp(mappingTable.ObjectID{i}, clearedMets));
        for j= 1:length(metIndices)
            metIndex= metIndices(j);
            if num2str(mappingTable.ChEBI(i)) ~= "NaN"
                model = addMIRIAMAnnotations(model, model.mets(metIndex),...
                    'ChEBI' , strcat('CHEBI:', num2str(mappingTable.ChEBI(i))));
            end
            if ~isempty(mappingTable.BIGG{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex),...
                    'bigg.metabolite' , mappingcTable.BIGG(i));
            end
            if num2str(mappingTable.ChemSpider(i)) ~= "NaN"
                model = addMIRIAMAnnotations(model, model.mets(metIndex),...
                    'ChemSpider', num2str(mappingTable.ChemSpider(i)));
            end
            if ~isempty(mappingTable.MetaNetX{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex), ...
                    'MetaNetX chemical' , mappingTable.MetaNetX(i));
            end
    
            if num2str(mappingTable.PubChem(i)) ~= "NaN"
                model = addMIRIAMAnnotations(model, model.mets(metIndex),...
                    'PubChem-compound', num2str(mappingTable.PubChem(i)));
            end
    
            if ~isempty(mappingTable.Seed{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex), ...
                    'SEED Compound' , mappingTable.Seed(i));
            end
    
            if ~isempty(mappingTable.InChI{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex), ...
                    'InChI' , mappingTable.InChI(i));
            end
            if ~isempty(mappingTable.InChI_Key{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex), ...
                    'InChIKey' , replace(mappingTable.InChI_Key(i), 'InChIKey=', ''));
            end
            if ~isempty(mappingTable.Kegg{i})
                model = addMIRIAMAnnotations(model, model.mets(metIndex), ...
                    'KEGG Compound' , mappingTable.Kegg(i));
            end
        end
    end
    
    model.metInChIString= model.metInChIString(1:length(model.mets));
    model.metKEGGID= model.metKEGGID(1:length(model.mets));
    model.metChEBIID= model.metChEBIID(1:length(model.mets));
    model.metPubChemID= model.metPubChemID(1:length(model.mets));
    
    for i= length(model.rxnECNumbers): length(model.rxns)
        model.rxnECNumbers(i)= {''};                        
    end


end
function [model, metHasMoved, originalMets]= moveRxnToComp(model, rxnIndex, compToMove, metHasMoved, originalMets)
    % this function moves a corresponding metabolites of a rxn to the
    % target compartment by manipulating the stoichiometric matrix
    % inputs:
    %   model:            original model
    %   rxnIndex:         index of the rxn in the model
    %   compToMove:       the target compartment that the rxn should be moved
    %   metHasMoved:      binary array that checkes weather the metabolite
    %                       has been moved to any compartment or not.
    %   original Mets     array of name of the metabolites without any compartment 
    % outputs:
    %   model:            updated model
    %   metHasMoved:      updateed metHasMoved
    %   originalMets:     updated originalMets

    % finding the metabolites participating in the rxn
    mets= find(model.S(:, rxnIndex));
    % loopeing through the mets of the rxn
    for i= 1:numel(mets)
        metIndex= mets(i);   
        compToMoveChar= model.comps{compToMove};
        % check weather the met has moved before
        % if it has not been moved yet, the changes are only applied to
        % itself
        if ~metHasMoved(metIndex)
            model.metComps(metIndex)= compToMove;
            model.mets(metIndex)= strcat(model.mets(metIndex), '[', compToMoveChar, ']');
            metHasMoved(metIndex)= true;
        % if the met has been moved before, the compartment is applied to
        % its duplicated version
        else
            if model.metComps(metIndex) ~= compToMove
                sCoef= model.S(metIndex, rxnIndex);
                newMet= strcat(originalMets(metIndex), '[', compToMoveChar, ']');

                % also other propertise of the metabolite should be included such as charge, etc.
                newMetName= model.metNames(metIndex); 
                % if the metabolite does not exist in the new compartment
                % it should be created
                if ~ismember(newMet, model.mets)
                    model= addMetabolite(model, newMet, 'metName', newMetName,...
                        'metFormula', model.metFormulas(metIndex),...
                        'ChEBIID', model.metChEBIID(metIndex), ...
                        'KEGGId', model.metKEGGID(metIndex), ...
                        'PubChemID', model.metPubChemID(metIndex), ...
                        'InChi', model.metInChIString(metIndex), ...
                        'Charge', model.metCharges(metIndex));
                    model.S(metIndex, rxnIndex)= 0;
                    model.S(end, rxnIndex)= sCoef;
                    model.metComps(end)= compToMove;
                    metHasMoved(end+1)= 1;
                    originalMets(end+1)= originalMets(metIndex);
                % if there is a duplicate version of the metabolite in the
                % target compartment, the changes are applied to it and
                % there is no more duplicating here.
                else
                    existingMetIndex= find(strcmp(newMet, model.mets));
                    model.S(metIndex, rxnIndex)= 0;
                    model.S(existingMetIndex, rxnIndex)= sCoef;
                end
            end
        end
    end
end
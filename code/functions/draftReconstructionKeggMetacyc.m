function draft= draftReconstructionKeggMetacyc(inputFastaFile)
    % reconstruction based on kegg and metacyc
    % inputs:
    %       inputFastaFile    the protein fasta file of the organism
    % output:
    %       draft             reconstructed draft

    keggModel=getKEGGModelForOrganism('chlorella_ohadii', inputFastaFile, ...
        'euk90_kegg105', 'output', ...
        true, false, false, false);
    metacycModel= getMetaCycModelForOrganism('chlorella_ohadii', inputFastaFile, ...
        false, false, false);
    draft=combineMetaCycKEGGModels(metacycModel,keggModel);
    draft.id= 'iCO';
    draft.name='Chlorella Ohadii';
end
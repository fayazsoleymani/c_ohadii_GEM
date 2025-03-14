initCobraToolbox;

%% model reconstruction from kegg and metacyc
proteinFastaFile= fullfile('..', 'data', 'ncbi', 'protein.faa');
draft= draftReconstructionKeggMetacyc(proteinFastaFile);

%%

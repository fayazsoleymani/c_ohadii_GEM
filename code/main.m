% model reconstruction from kegg and metacyc
proteinFastaFile= fullfile('..', 'data', 'ncbi', 'protein.faa');
draft= draftReconstructionKeggMetacyc(proteinFastaFile);

% exchange reactions
exchange2lbDict= dictionary({'WATER'}, -10, {'PROTON'}, -10, {'Pi'}, -10, ...
    {'AMMONIUM'}, -1, {'SULFATE'}, -10, {'FE+2'}, -10, {'CARBON-DIOXIDE'}, ...
    -2, {'MG+2'}, -10, {'C00205'}, -80, {'OXYGEN-MOLECULE'}, -10);
model= modifyExchangeRxns(model, exchange2lbDict);

model= addBiomassRxn(model, fullfile('..', 'data', 'biomass_all.tsv'));


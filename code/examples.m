%% model reconstruction from kegg and metacyc
proteinFastaFile= fullfile('..', 'data', 'ncbi', 'protein.faa');
draft= draftReconstructionKeggMetacyc(proteinFastaFile);

%% adding exchange reactions
exchange2lbDict= dictionary({'WATER'}, -10, {'PROTON'}, -10, {'Pi'}, -10, ...
    {'AMMONIUM'}, -1, {'SULFATE'}, -10, {'FE+2'}, -10, {'CARBON-DIOXIDE'}, ...
    -2, {'MG+2'}, -10, {'C00205'}, -80, {'OXYGEN-MOLECULE'}, -10, ...
    {'ACET'}, -10);
model= modifyExchangeRxns(model, exchange2lbDict);

%% adding biomass reactions to the model
model= addBiomassRxn(model, fullfile('..', 'data', 'biomass_all.tsv'));

%% modifying the GPR rules for the lipid module
blastResultsPath= '../data/genes/blast_at_cohadii.csv';
model= modifyGPRsBlast(model, blastResultsPath);

%% compartmentalization first part
solutionBefore= checkBiomassProduction(model, 'one', false);

GSSfilePath= fullfile('..', 'data', 'genes', 'GSS_file_new.csv');
GSS= constructGSS(GSSfilePath);

predefinedCompartments= readtable(fullfile('..', 'data', 'genes', ...
    'predefined_rxn_comps_clean.csv'), 'Delimiter', ',', ...
    'VariableNamingRule', 'preserve');

model= compartmentalize_first(model, GSS, 0.25, 'one', predefinedCompartments);

solutionAfter= checkBiomassProduction(model, 'one', false);

%%


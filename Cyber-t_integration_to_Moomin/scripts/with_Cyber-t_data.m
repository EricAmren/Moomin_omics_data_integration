clear;

% A MILP-solver needs to be set for COBRA TB
changeCobraSolver('ibm_cplex','milp');

% WITH CYBER-T DATA

% load the model and the expression data
model = readSBML('example/iJO1366.xml',10);
%model = readSBML('example/Ec_core_flux1.xml',10);

cyber_t_expression = tdfread('/home/ericcumunel/Documents/M2_Internship/microarray_VS_RNAseq/Moomin_input/cyber_t_res.tsv')

% solve for one solution using stoichiometric contraints
[cyber_t_model_stoich, cyber_t_MILPsolution_stoich, cyber_t_MILPproblem_stoich] = moomin(model, cyber_t_expression);


% solve for up to 10 solutions with only topological contraints
[cyber_t_model_topo, cyber_t_MILPsolution_topo, cyber_t_MILPproblem_topo] = moomin(model, cyber_t_expression, 'stoichiometry', 0, 'enumerate', 10);

% write the standard output file
writeOutput(cyber_t_model_stoich, 'example/cyber_t_output_stoich.json', 'type', 'json');


% write the standard output for input colours
writeOutput(cyber_t_model_topo, 'example/cyber_t_input_topo.json','type','json');

writeSBML(model,'result_cyber_t_model.xml')
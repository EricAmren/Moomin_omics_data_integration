clear;

% A MILP-solver needs to be set for COBRA TB
changeCobraSolver('ibm_cplex','milp');

% WITH CYBER-T DATA

model_for_cyber_t = readSBML('models/iJO1366.xml',10);
cyber_t_expression = tdfread('Moomin_input/cyber_t_res.tsv');
[cyber_t_model_stoich, cyber_t_MILPsolution_stoich, cyber_t_MILPproblem_stoich] = moomin(model_for_cyber_t, cyber_t_expression);
[cyber_t_model_topo, cyber_t_MILPsolution_topo, cyber_t_MILPproblem_topo] = moomin(model_for_cyber_t, cyber_t_expression, 'stoichiometry', 0, 'enumerate', 10);
writeOutput(cyber_t_model_stoich, 'output/cyber_t_output_stoich.json', 'type', 'json');
writeOutput(cyber_t_model_topo, 'output/cyber_t_input_topo.json','type','json');
writeSBML(model_for_cyber_t,'output/result_cyber_t_model.xml')

% WITH EBSEQ DATA
model_for_EBseq = readSBML('models/iJO1366.xml',10);
EBseq_expression = tdfread('Moomin_input/EBseq_res.tsv');
[EBseq_model_stoich, EBseq_MILPsolution_stoich, EBseq_MILPproblem_stoich] = moomin(model_for_EBseq, EBseq_expression);
[EBseq_model_topo, EBseq_MILPsolution_topo, EBseq_MILPproblem_topo] = moomin(model_for_EBseq, EBseq_expression, 'stoichiometry', 0, 'enumerate', 10);
writeOutput(EBseq_model_stoich, 'output/EBseq_output_stoich.json', 'type', 'json');
writeOutput(EBseq_model_topo, 'output/EBseq_input_topo.json','type','json');
writeSBML(model_for_EBseq,'output/result_EBseq_model.xml')

% Comparison

concordance_matrix = cyber_t_model_stoich.outputColours==EBseq_model_stoich.outputColours;
d = num2cell(c), cyber_t_model_stoich
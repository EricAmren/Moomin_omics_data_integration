clear;

% A MILP-solver needs to be set for COBRA TB
changeCobraSolver('ibm_cplex','milp');

model = readSBML('example/iJO1366.xml',10);

EBseq_expression = tdfread('/home/ericcumunel/Documents/M2_Internship/microarray_VS_RNAseq/Moomin_input/EBseq_res.tsv');
[EBseq_model_stoich, EBseq_MILPsolution_stoich, EBseq_MILPproblem_stoich] = moomin(model, EBseq_expression);
[EBseq_model_topo, EBseq_MILPsolution_topo, EBseq_MILPproblem_topo] = moomin(model, EBseq_expression, 'stoichiometry', 0, 'enumerate', 10);
writeOutput(EBseq_model_stoich, 'example/EBseq_output_stoich.json', 'type', 'json');
writeOutput(EBseq_model_topo, 'example/EBseq_input_topo.json','type','json');

writeSBML(model,'result_EBseq_model.xml')
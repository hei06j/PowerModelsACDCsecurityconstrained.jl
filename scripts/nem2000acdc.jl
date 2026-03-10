using Pkg
Pkg.activate("./")

using Ipopt
using JuMP
using PowerModels
using PowerModelsACDC
using PowerModelsSecurityConstrained
using PowerModelsACDCsecurityconstrained


const _PM = PowerModels
const _PMACDC = PowerModelsACDC
const _PMSC = PowerModelsSecurityConstrained
const _PMSCACDC = PowerModelsACDCsecurityconstrained


nlp_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0) 

file = "./test/data/matpower/snem2000_acdc.m"
data = parse_file(file)
_PMSCACDC.fix_scopf_data_issues!(data)
_PMACDC.process_additional_data!(data)
setting = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => true) 

# test for only 1 branch contingency
data["gen_contingencies"] = [data["gen_contingencies"][1]]
# data["gen_contingencies"] = data["gen_contingencies"][1:10]
data["branch_contingencies"] = data["branch_contingencies"][1:10]
data["branchdc_contingencies"] = []
data["convdc_contingencies"] = []

_PMSCACDC.silence()

result = _PMSCACDC.run_scopf_acdc_contingencies(data, _PM.ACPPowerModel, _PM.ACPPowerModel, _PMSCACDC.run_scopf, nlp_solver, nlp_solver, setting)


# network = deepcopy(data)
# model_type_scopf = _PM.ACPPowerModel
# model_type_filter = _PM.ACPPowerModel
# run_scopf_prob = _PMSCACDC.run_scopf
# optimizer_scopf = nlp_solver
# optimizer_filter = nlp_solver
# max_iter = 100
# time_limit = Inf
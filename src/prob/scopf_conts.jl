"""
This formulation is best used in conjunction with the contingency filters that find
violated contingencies in integrated ACDC grid.

"""

function run_scopf(data, model_constructor, solver; kwargs...)
    # _PMACDC.process_additional_data!(data)
    return _PM.solve_model(data, model_constructor, solver, build_scopf; ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PMACDC.add_ref_dcgrid!, _PMACDC.ref_add_pst!, _PMACDC.ref_add_sssc!,_PMSC.ref_c1!], multinetwork=true, kwargs...)
end


# no data, so no further templating is needed, constraint goes directly to the formulations
function constraint_power_balance_ac(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    bus = _PM.ref(pm, nw, :bus, i)
    bus_arcs = _PM.ref(pm, nw, :bus_arcs, i)
    bus_arcs_pst = _PM.ref(pm, nw, :bus_arcs_pst, i)
    bus_arcs_sssc = _PM.ref(pm, nw, :bus_arcs_sssc, i)
    bus_arcs_sw = _PM.ref(pm, nw, :bus_arcs_sw, i)
    bus_gens = _PM.ref(pm, nw, :bus_gens, i)
    bus_loads = _PM.ref(pm, nw, :bus_loads, i)
    bus_shunts = _PM.ref(pm, nw, :bus_shunts, i)
    bus_storage = _PM.ref(pm, nw, :bus_storage, i)
    bus_convs_ac = _PM.ref(pm, nw, :bus_convs_ac, i)

    bus_pd = Dict(k => _PM.ref(pm, nw, :load, k, "pd") for k in bus_loads)
    bus_qd = Dict(k => _PM.ref(pm, nw, :load, k, "qd") for k in bus_loads)

    bus_gs = Dict(k => _PM.ref(pm, nw, :shunt, k, "gs") for k in bus_shunts)
    bus_bs = Dict(k => _PM.ref(pm, nw, :shunt, k, "bs") for k in bus_shunts)

    constraint_power_balance_ac(pm, nw, i, bus_arcs, bus_arcs_pst, bus_arcs_sssc, bus_convs_ac, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
end

function constraint_power_balance_ac(pm::_PM.AbstractACPModel, n::Int, i::Int, bus_arcs, bus_arcs_pst, bus_arcs_sssc, bus_convs_ac, bus_arcs_sw, bus_gens, bus_storage, bus_loads, bus_gs, bus_bs)
    vm   = _PM.var(pm, n, :vm, i)
    p    = _PM.var(pm, n, :p)
    q    = _PM.var(pm, n, :q)
    # ppst    = _PM.var(pm, n,    :ppst)
    # qpst    = _PM.var(pm, n,    :qpst)
    pg   = _PM.var(pm, n,   :pg)
    qg   = _PM.var(pm, n,   :qg)
    pconv_grid_ac = _PM.var(pm, n,  :pconv_tf_fr)
    qconv_grid_ac = _PM.var(pm, n,  :qconv_tf_fr)
    # pflex = _PM.var(pm, n, :pflex)
    # qflex = _PM.var(pm, n, :qflex)
    # ps    = _PM.var(pm, n, :ps)
    # qs    = _PM.var(pm, n, :qs)
    # psssc    = _PM.var(pm, n,    :psssc)
    # qsssc    = _PM.var(pm, n,    :qsssc)


    cstr_p = JuMP.@constraint(pm.model,
        sum(p[a] for a in bus_arcs)
        # + sum(ppst[a] for a in bus_arcs_pst)
        # + sum(psssc[a] for a in bus_arcs_sssc) 
        + sum(pconv_grid_ac[c] for c in bus_convs_ac)
        ==
        sum(pg[g] for g in bus_gens)
        # - sum(ps[s] for s in bus_storage)
        # - sum(pflex[d] for d in bus_loads)
        - sum(gs for (i,gs) in bus_gs)*vm^2
    )
    cstr_q = JuMP.@constraint(pm.model,
        sum(q[a] for a in bus_arcs)
        # + sum(qpst[a] for a in bus_arcs_pst)
        # + sum(qsssc[a] for a in bus_arcs_sssc) 
        + sum(qconv_grid_ac[c] for c in bus_convs_ac)
        ==
        sum(qg[g] for g in bus_gens)
        # - sum(qs[s] for s in bus_storage)
        # - sum(qflex[d] for d in bus_loads)
        + sum(bs for (i,bs) in bus_bs)*vm^2
    )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :bus, i)[:lam_kcl_r] = cstr_p
        _PM.sol(pm, n, :bus, i)[:lam_kcl_i] = cstr_q
    end
end

function constraint_power_balance_dc(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    bus_arcs_dcgrid = _PM.ref(pm, nw, :bus_arcs_dcgrid, i)
    bus_convs_dc = _PM.ref(pm, nw, :bus_convs_dc, i)
    bus_gens_dc = [] #_PM.ref(pm, nw, :bus_gens_dc, i)
    pd = _PM.ref(pm, nw, :busdc, i)["Pdc"]
    constraint_power_balance_dc(pm, nw, i, bus_arcs_dcgrid, bus_convs_dc, bus_gens_dc, pd)
end

function constraint_power_balance_dc(pm::_PM.AbstractPowerModel, n::Int, i::Int, bus_arcs_dcgrid, bus_convs_dc, bus_gens_dc, pd)
    p_dcgrid = _PM.var(pm, n, :p_dcgrid)
    pconv_dc = _PM.var(pm, n, :pconv_dc)
    # pgdc = _PM.var(pm, n, :pgdc)

    cstr_p = JuMP.@constraint(pm.model, 
            sum(p_dcgrid[a] for a in bus_arcs_dcgrid) + 
            sum(pconv_dc[c] for c in bus_convs_dc) 
            == (-pd) 
            # + sum(pgdc[g] for g in bus_gens_dc)
            )

    if _IM.report_duals(pm)
        _PM.sol(pm, n, :busdc, i)[:lam_kcl_r] = cstr_p
    end
end

function constraint_converter_current(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    Vmax = conv["Vmmax"]
    Imax = conv["Imax"]
    constraint_converter_current(pm, nw, i, Vmax, Imax)
end

function constraint_converter_current(pm::_PM.AbstractACPModel, n::Int, i::Int, Umax, Imax)
    vmc = _PM.var(pm, n, :vmc, i)
    pconv_ac = _PM.var(pm, n, :pconv_ac, i)
    qconv_ac = _PM.var(pm, n, :qconv_ac, i)
    iconv = _PM.var(pm, n, :iconv_ac, i)

    JuMP.@NLconstraint(pm.model, pconv_ac^2 + qconv_ac^2 == vmc^2 * iconv^2)
end

function constraint_conv_transformer(pm::_PM.AbstractPowerModel, i::Int; nw::Int=_PM.nw_id_default)
    conv = _PM.ref(pm, nw, :convdc, i)
    constraint_conv_transformer(pm, nw, i, conv["rtf"], conv["xtf"], conv["busac_i"], conv["tm"], Bool(conv["transformer"]))
end


function constraint_conv_transformer(pm::_PM.AbstractACPModel, n::Int, i::Int, rtf, xtf, acbus, tm, transformer)
    ptf_fr = _PM.var(pm, n, :pconv_tf_fr, i)
    qtf_fr = _PM.var(pm, n, :qconv_tf_fr, i)
    ptf_to = _PM.var(pm, n, :pconv_tf_to, i)
    qtf_to = _PM.var(pm, n, :qconv_tf_to, i)

    vm = _PM.var(pm, n, :vm, acbus)
    va = _PM.var(pm, n, :va, acbus)
    vmf = _PM.var(pm, n, :vmf, i)
    vaf = _PM.var(pm, n, :vaf, i)

    ztf = rtf + im*xtf
    if transformer
        ytf = 1/(rtf + im*xtf)
        gtf = real(ytf)
        btf = imag(ytf)
        gtf_sh = 0
        c1, c2, c3, c4 = ac_power_flow_constraints(pm.model, gtf, btf, gtf_sh, vm, vmf, va, vaf, ptf_fr, ptf_to, qtf_fr, qtf_to, tm)
    else
        JuMP.@constraint(pm.model, ptf_fr + ptf_to == 0)
        JuMP.@constraint(pm.model, qtf_fr + qtf_to == 0)
        JuMP.@constraint(pm.model, va == vaf)
        JuMP.@constraint(pm.model, vm == vmf)
    end
end

function ac_power_flow_constraints(model, g, b, gsh_fr, vm_fr, vm_to, va_fr, va_to, p_fr, p_to, q_fr, q_to, tm)
    c1 = JuMP.@NLconstraint(model, p_fr ==  g/(tm^2)*vm_fr^2 + -g/(tm)*vm_fr*vm_to * cos(va_fr-va_to) + -b/(tm)*vm_fr*vm_to*sin(va_fr-va_to))
    c2 = JuMP.@NLconstraint(model, q_fr == -b/(tm^2)*vm_fr^2 +  b/(tm)*vm_fr*vm_to * cos(va_fr-va_to) + -g/(tm)*vm_fr*vm_to*sin(va_fr-va_to))
    c3 = JuMP.@NLconstraint(model, p_to ==  g*vm_to^2 + -g/(tm)*vm_to*vm_fr  *    cos(va_to - va_fr)     + -b/(tm)*vm_to*vm_fr    *sin(va_to - va_fr))
    c4 = JuMP.@NLconstraint(model, q_to == -b*vm_to^2 +  b/(tm)*vm_to*vm_fr  *    cos(va_to - va_fr)     + -g/(tm)*vm_to*vm_fr    *sin(va_to - va_fr))
    return c1, c2, c3, c4
end


function build_scopf(pm::_PM.AbstractPowerModel)
    
    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)
    _PM.variable_branch_transform(pm, nw=0)
    _PM.constraint_model_voltage(pm, nw=0)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=0)

    _PMACDC.variable_active_dcbranch_flow(pm, nw=0)
    _PMACDC.variable_dcbranch_current(pm, nw=0)
    _PMACDC.variable_dc_converter(pm, nw=0)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=0)
    _PMACDC.constraint_voltage_dc(pm, nw=0)
    
    
    for i in _PM.ids(pm, nw=0, :ref_buses)
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :bus)
        # _PMACDC.constraint_power_balance_ac(pm, i, nw=0)
        constraint_power_balance_ac(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branch)
        constraint_ohms_y_oltc_pst_from(pm, i, nw=0)
        constraint_ohms_y_oltc_pst_to(pm, i, nw=0)
        _PM.constraint_voltage_angle_difference(pm, i, nw=0)
        _PM.constraint_thermal_limit_from(pm, i, nw=0)
        _PM.constraint_thermal_limit_to(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :busdc)                          
        # _PMACDC.constraint_power_balance_dc(pm, i, nw=0)
        constraint_power_balance_dc(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :convdc)
        _PMACDC.constraint_converter_losses(pm, i, nw=0)
        # _PMACDC.constraint_converter_current(pm, i, nw=0)
        # _PMACDC.constraint_conv_transformer(pm, i, nw=0)
        constraint_converter_current(pm, i, nw=0)
        constraint_conv_transformer(pm, i, nw=0)
        _PMACDC.constraint_conv_reactor(pm, i, nw=0)
        _PMACDC.constraint_conv_filter(pm, i, nw=0)
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1
            _PMACDC.constraint_conv_firing_angle(pm, i, nw=0)
        end
        constraint_pvdc_droop_control_linear(pm, i, nw=0)
    end
    
    contigency_ids = [id for id in _PM.nw_ids(pm) if id != 0]
    
    for nw in contigency_ids
        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)    
        _PM.variable_gen_power(pm, nw=nw, bounded=false)       
        _PM.variable_branch_power(pm, nw=nw)
        _PM.variable_branch_transform(pm, nw=nw)
        _PM.constraint_model_voltage(pm, nw=nw)

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=nw)

        _PMACDC.variable_active_dcbranch_flow(pm, nw=nw)
        _PMACDC.variable_dcbranch_current(pm, nw=nw)
        _PMACDC.variable_dc_converter(pm, nw=nw)
        _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw)
        _PMACDC.constraint_voltage_dc(pm, nw=nw)

        _PMSC.variable_c1_response_delta(pm, nw=nw)

        for i in _PM.ids(pm, nw=nw, :ref_buses)
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = _PM.ref(pm, :gen_buses, nw=nw)
        for i in _PM.ids(pm, :bus, nw=nw)
            _PMACDC.constraint_power_balance_ac(pm, i, nw=nw)

            # if a bus has active generators, fix the voltage magnitude to the base case
            if i in gen_buses
                _PMSC.constraint_c1_voltage_magnitude_link(pm, i, nw_1=0, nw_2=nw)
            end
        end

        response_gens = _PM.ref(pm, :response_gens, nw=nw)
        for (i,gen) in _PM.ref(pm, :gen, nw=nw)
            pg_base = _PM.var(pm, :pg, i, nw=0)

            # setup the linear response function or fix value to base case
            if i in response_gens
                _PMSC.constraint_c1_gen_power_real_response(pm, i, nw_1=0, nw_2=nw)
            else
                _PMSC.constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end

        for i in _PM.ids(pm,  nw=nw, :branch)
            constraint_ohms_y_oltc_pst_from(pm, i, nw=nw)
            constraint_ohms_y_oltc_pst_to(pm, i, nw=nw)
            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)
            _PM.constraint_thermal_limit_from(pm, i, nw=nw)
            _PM.constraint_thermal_limit_to(pm, i, nw=nw)
        end

        for i in _PM.ids(pm, nw=nw, :busdc)                        
            _PMACDC.constraint_power_balance_dc(pm, i, nw=nw)
        end                                          
        for i in _PM.ids(pm, nw=nw, :branchdc)
            _PMACDC.constraint_ohms_dc_branch(pm, i, nw=nw)
        end                                          
        for i in _PM.ids(pm, nw=nw, :convdc)          
            _PMACDC.constraint_converter_losses(pm, i, nw=nw)
            _PMACDC.constraint_converter_current(pm, i, nw=nw)
            _PMACDC.constraint_conv_transformer(pm, i, nw=nw)
            _PMACDC.constraint_conv_reactor(pm, i, nw=nw)
            _PMACDC.constraint_conv_filter(pm, i, nw=nw)
            if pm.ref[:it][:pm][:nw][nw][:convdc][i]["islcc"] == 1
                _PMACDC.constraint_conv_firing_angle(pm, i, nw=nw)
            end
            constraint_pvdc_droop_control_linear(pm, i, nw=nw)                                      
        end
 
    end

    objective_min_fuel_cost_scopf(pm)
    
end



"""
An SCOPF multi-period soft formulation for integrated HVAC and HVDC grid. It includes
slack variables for AC and DC grid power balance and line thermal limit constraints, 
which are minimized in the objective function.

This formulation is best used in conjunction with the contingency filters that find
violated contingencies in integrated HVAC and HVDC grid.

"""

function run_scopf_soft(data, model_constructor, solver; kwargs...)
    # _PMACDC.process_additional_data!(data)
    return _PM.solve_model(data, model_constructor, solver, build_scopf_soft; ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PMACDC.add_ref_dcgrid!, _PMSC.ref_c1!], multinetwork=true, kwargs...) 
end

function build_scopf_soft(pm::_PM.AbstractPowerModel)
    
    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)
    _PM.variable_branch_transform(pm, nw=0)
    _PM.constraint_model_voltage(pm, nw=0)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=0)

    _PMACDC.variable_active_dcbranch_flow(pm, nw=0)       
    _PMACDC.variable_dcbranch_current(pm, nw=0)
    _PMACDC.variable_dc_converter(pm, nw=0)                        
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=0) 
    _PMACDC.constraint_voltage_dc(pm, nw=0)    
   
    variables_slacks(pm, nw=0)

    for i in _PM.ids(pm, nw=0, :ref_buses)                  
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :bus)                        
        constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branch)                     
        constraint_ohms_y_oltc_pst_from(pm, i, nw=0)
        constraint_ohms_y_oltc_pst_to(pm, i, nw=0)
        _PM.constraint_voltage_angle_difference(pm, i, nw=0)
        constraint_thermal_limit_from_soft(pm, i, nw=0)
        constraint_thermal_limit_to_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :busdc)                                        
        constraint_power_balance_dc_soft(pm, i, nw=0)                                      
    end

    for i in _PM.ids(pm, nw=0, :branchdc)                                             
        constraint_ohms_dc_branch_soft(pm, i, nw=0)                                 
    end 

    for i in _PM.ids(pm, nw=0, :convdc)                                                
        _PMACDC.constraint_converter_losses(pm, i, nw=0)                                                            
        _PMACDC.constraint_conv_transformer(pm, i, nw=0)                               
        _PMACDC.constraint_conv_reactor(pm, i, nw=0)                                 
        _PMACDC.constraint_conv_filter(pm, i, nw=0)
        constraint_converter_current(pm, i, nw=0)                                     
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1                
            _PMACDC.constraint_conv_firing_angle(pm, i, nw=0)                                     
        end
        constraint_pvdc_droop_control_linear(pm, i, nw=0)                                                                        
    end                                                 

    contigency_ids = [id for id in _PM.nw_ids(pm) if id != 0]         
    
    for nw in contigency_ids

        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)
        _PM.variable_gen_power(pm, nw=nw, bounded=false)
        _PM.variable_branch_power(pm, nw=nw)
        _PM.variable_branch_transform(pm, nw=nw)
        _PM.constraint_model_voltage(pm, nw=nw)

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=nw)
        _PMSC.variable_c1_response_delta(pm, nw=nw)

        _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=false)       
        _PMACDC.variable_dcbranch_current(pm, nw=nw)
        _PMACDC.variable_dc_converter(pm, nw=nw)                        
        _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw)
        _PMACDC.constraint_voltage_dc(pm, nw=nw)

        variables_slacks(pm, nw=nw)

        for i in _PM.ids(pm, nw=nw, :ref_buses)           
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = _PM.ref(pm, nw=nw, :gen_buses)
        for i in _PM.ids(pm, nw=nw, :bus)                 
            constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=nw)       

            # if a bus has active generators, fix the voltage magnitude to the base case
            if i in gen_buses
                _PMSC.constraint_c1_voltage_magnitude_link(pm, i, nw_1=0, nw_2=nw)
            end
        end

        response_gens = _PM.ref(pm, :response_gens, nw=nw)
        for (i,gen) in _PM.ref(pm, :gen, nw=nw)
            pg_base = _PM.var(pm, :pg, i, nw=0)

            # setup the linear response function or fix value to base case
            if i in response_gens
                _PMSC.constraint_c1_gen_power_real_response(pm, i, nw_1=0, nw_2=nw)
            else
                _PMSC.constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end
  
        for i in _PM.ids(pm, nw=nw, :branch)                     
            constraint_ohms_y_oltc_pst_from(pm, i, nw=nw)
            constraint_ohms_y_oltc_pst_to(pm, i, nw=nw)
            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)
            constraint_thermal_limit_from_soft(pm, i, nw=nw)
            constraint_thermal_limit_to_soft(pm, i, nw=nw)
        end

        for i in _PM.ids(pm, nw=nw, :busdc)                                               
            constraint_power_balance_dc_soft(pm, i, nw=nw)                                     
        end                                                                               
        for i in _PM.ids(pm, nw=nw, :branchdc)                                             
            constraint_ohms_dc_branch_soft(pm, i, nw=nw)                                
        end                                                                                
        for i in _PM.ids(pm, nw=nw, :convdc)                                               
            _PMACDC.constraint_converter_losses(pm, i, nw=nw)                              
            constraint_converter_current(pm, i, nw=nw)                             
            _PMACDC.constraint_conv_transformer(pm, i, nw=nw)                              
            _PMACDC.constraint_conv_reactor(pm, i, nw=nw)                                  
            _PMACDC.constraint_conv_filter(pm, i, nw=nw)                                   
            if pm.ref[:it][:pm][:nw][nw][:convdc][i]["islcc"] == 1            
                _PMACDC.constraint_conv_firing_angle(pm, i, nw=nw)                                    
            end
            constraint_pvdc_droop_control_linear(pm, i, nw=nw)                                                                            
        end
 
    end
    objective_min_fuel_cost_scopf_soft(pm)
       
end

function run_scopf_soft_frq(data, model_constructor, solver; kwargs...)
    # _PMACDC.process_additional_data!(data)
    return _PM.solve_model(data, model_constructor, solver, build_scopf_soft_frq; ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PMACDC.add_ref_dcgrid!, _PMSC.ref_c1!], multinetwork=true, kwargs...) 
end

function build_scopf_soft_frq(pm::_PM.AbstractPowerModel)
    
    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)
    _PM.variable_branch_transform(pm, nw=0)
    _PM.constraint_model_voltage(pm, nw=0)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=0)

    _PMACDC.variable_active_dcbranch_flow(pm, nw=0)       
    _PMACDC.variable_dcbranch_current(pm, nw=0)
    _PMACDC.variable_dc_converter(pm, nw=0)                        
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=0) 
    _PMACDC.constraint_voltage_dc(pm, nw=0)    
   
    variables_slacks(pm, nw=0)

    # variable_area_frequency(pm, nw=0)

    for i in _PM.ids(pm, nw=0, :ref_buses)                  
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :bus)                        
        constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branch)                     
        constraint_ohms_y_oltc_pst_from(pm, i, nw=0)
        constraint_ohms_y_oltc_pst_to(pm, i, nw=0)
        _PM.constraint_voltage_angle_difference(pm, i, nw=0)
        constraint_thermal_limit_from_soft(pm, i, nw=0)
        constraint_thermal_limit_to_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :busdc)                                        
        constraint_power_balance_dc_soft(pm, i, nw=0)                                      
    end

    for i in _PM.ids(pm, nw=0, :branchdc)                                             
        constraint_ohms_dc_branch_soft(pm, i, nw=0)                                 
    end 

    for i in _PM.ids(pm, nw=0, :convdc)                                                
        _PMACDC.constraint_converter_losses(pm, i, nw=0)                                                            
        _PMACDC.constraint_conv_transformer(pm, i, nw=0)                               
        _PMACDC.constraint_conv_reactor(pm, i, nw=0)                                 
        _PMACDC.constraint_conv_filter(pm, i, nw=0)
        constraint_converter_current(pm, i, nw=0)                                     
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1                
            _PMACDC.constraint_conv_firing_angle(pm, i, nw=0)                                     
        end
        constraint_pvdc_droop_control_linear(pm, i, nw=0)                                                                        
    end                                                 

    contigency_ids = [id for id in _PM.nw_ids(pm) if id != 0]         
    
    for nw in contigency_ids

        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)
        _PM.variable_gen_power(pm, nw=nw, bounded=false)
        _PM.variable_branch_power(pm, nw=nw)
        _PM.variable_branch_transform(pm, nw=nw)
        _PM.constraint_model_voltage(pm, nw=nw)

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=nw)
        _PMSC.variable_c1_response_delta(pm, nw=nw)

        _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=false)       
        _PMACDC.variable_dcbranch_current(pm, nw=nw)
        _PMACDC.variable_dc_converter(pm, nw=nw)                        
        _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw)
        _PMACDC.constraint_voltage_dc(pm, nw=nw)

        variables_slacks(pm, nw=nw)

        # variable_area_frequency(pm, nw=nw)

        for i in _PM.ids(pm, nw=nw, :ref_buses)           
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = _PM.ref(pm, nw=nw, :gen_buses)
        for i in _PM.ids(pm, nw=nw, :bus)                 
            constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=nw)       

            # if a bus has active generators, fix the voltage magnitude to the base case
            if i in gen_buses
                _PMSC.constraint_c1_voltage_magnitude_link(pm, i, nw_1=0, nw_2=nw)
            end
        end

        response_gens = _PM.ref(pm, :response_gens, nw=nw)
        for (i,gen) in _PM.ref(pm, :gen, nw=nw)
            pg_base = _PM.var(pm, :pg, i, nw=0)

            # setup the linear response function or fix value to base case
            if i in response_gens
                _PMSC.constraint_c1_gen_power_real_response(pm, i, nw_1=0, nw_2=nw)
            else
                _PMSC.constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end
  
        for i in _PM.ids(pm, nw=nw, :branch)                     
            constraint_ohms_y_oltc_pst_from(pm, i, nw=nw)
            constraint_ohms_y_oltc_pst_to(pm, i, nw=nw)
            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)
            constraint_thermal_limit_from_soft(pm, i, nw=nw)
            constraint_thermal_limit_to_soft(pm, i, nw=nw)
        end

        for i in _PM.ids(pm, nw=nw, :busdc)                                               
            constraint_power_balance_dc_soft(pm, i, nw=nw)                                     
        end                                                                               
        for i in _PM.ids(pm, nw=nw, :branchdc)                                             
            constraint_ohms_dc_branch_soft(pm, i, nw=nw)                                
        end                                                                                
        for i in _PM.ids(pm, nw=nw, :convdc)                                               
            _PMACDC.constraint_converter_losses(pm, i, nw=nw)                              
            constraint_converter_current(pm, i, nw=nw)                             
            _PMACDC.constraint_conv_transformer(pm, i, nw=nw)                              
            _PMACDC.constraint_conv_reactor(pm, i, nw=nw)                                  
            _PMACDC.constraint_conv_filter(pm, i, nw=nw)                                   
            if pm.ref[:it][:pm][:nw][nw][:convdc][i]["islcc"] == 1            
                _PMACDC.constraint_conv_firing_angle(pm, i, nw=nw)                                    
            end
            constraint_pvdc_droop_control_linear(pm, i, nw=nw)                                                                            
        end
 
    end
    objective_min_fuel_cost_scopf_soft(pm)
       
end


"""
An SCOPF multi-period soft formulation for integrated HVAC and HVDC grid. It includes
slack variables for AC and DC grid power balance and line thermal limit constraints, 
which are minimized in the objective function.

This formulation is best used in conjunction with the contingency filters that find
violated contingencies in integrated HVAC and HVDC grid.

"""

function run_scopf_soft_smooth(data, model_constructor, solver; kwargs...)
    # _PMACDC.process_additional_data!(data)
    return _PM.solve_model(data, model_constructor, solver, build_scopf_soft_smooth; ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PMACDC.add_ref_dcgrid!, _PMSC.ref_c1!], multinetwork=true, kwargs...) 
end

function build_scopf_soft_smooth(pm::_PM.AbstractPowerModel)
    
    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)
    _PM.variable_branch_transform(pm, nw=0)
    _PM.constraint_model_voltage(pm, nw=0)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=0)

    _PMACDC.variable_active_dcbranch_flow(pm, nw=0)       
    _PMACDC.variable_dcbranch_current(pm, nw=0)
    _PMACDC.variable_dc_converter(pm, nw=0)                        
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=0) 
    _PMACDC.constraint_voltage_dc(pm, nw=0)    
   
    variables_slacks(pm, nw=0)

    for i in _PM.ids(pm, nw=0, :ref_buses)                  
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :bus)                        
        constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branch)                     
        constraint_ohms_y_oltc_pst_from(pm, i, nw=0)
        constraint_ohms_y_oltc_pst_to(pm, i, nw=0)
        _PM.constraint_voltage_angle_difference(pm, i, nw=0)
        constraint_thermal_limit_from_soft(pm, i, nw=0)
        constraint_thermal_limit_to_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :busdc)                                        
        constraint_power_balance_dc_soft(pm, i, nw=0)                                      
    end

    for i in _PM.ids(pm, nw=0, :branchdc)                                             
        constraint_ohms_dc_branch_soft(pm, i, nw=0)                                 
    end 

    for i in _PM.ids(pm, nw=0, :convdc)                                                
        _PMACDC.constraint_converter_losses(pm, i, nw=0)                                                            
        _PMACDC.constraint_conv_transformer(pm, i, nw=0)                               
        _PMACDC.constraint_conv_reactor(pm, i, nw=0)                                 
        _PMACDC.constraint_conv_filter(pm, i, nw=0)
        constraint_converter_current(pm, i, nw=0)                                     
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1                
            _PMACDC.constraint_conv_firing_angle(pm, i, nw=0)                                     
        end
        constraint_pvdc_droop_control_smooth(pm, i, nw=0)                                                                        
    end                                                 

    contigency_ids = [id for id in _PM.nw_ids(pm) if id != 0]         
    
    for nw in contigency_ids

        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)
        _PM.variable_gen_power(pm, nw=nw, bounded=false)
        _PM.variable_branch_power(pm, nw=nw)
        _PM.variable_branch_transform(pm, nw=nw)
        _PM.constraint_model_voltage(pm, nw=nw)

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=nw)
        _PMSC.variable_c1_response_delta(pm, nw=nw)

        _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=false)       
        _PMACDC.variable_dcbranch_current(pm, nw=nw)
        _PMACDC.variable_dc_converter(pm, nw=nw)                        
        _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw)
        _PMACDC.constraint_voltage_dc(pm, nw=nw)
          
        variables_slacks(pm, nw=nw)

        for i in _PM.ids(pm, nw=nw, :ref_buses)           
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = _PM.ref(pm, nw=nw, :gen_buses)
        for i in _PM.ids(pm, nw=nw, :bus)                 
            constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=nw)       

            # for a bus with active generators, fix voltage magnitude to base case until limits are reached, then switch PV/PQ bus
            if i in gen_buses
                constraint_gen_power_reactive_response_smooth(pm, i, nw_1=0, nw_2=nw)
            end
        end

        response_gens = _PM.ref(pm, :response_gens, nw=nw)
        for (i,gen) in _PM.ref(pm, :gen, nw=nw)
            pg_base = _PM.var(pm, :pg, i, nw=0)

            # setup the smooth response function until limits are reached or fix value to base case
            if i in response_gens
                constraint_gen_power_real_response_smooth(pm, i, nw_1=0, nw_2=nw)
            else
                _PMSC.constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end
  
        for i in _PM.ids(pm, nw=nw, :branch)                     
            constraint_ohms_y_oltc_pst_from(pm, i, nw=nw)
            constraint_ohms_y_oltc_pst_to(pm, i, nw=nw)
            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)
            constraint_thermal_limit_from_soft(pm, i, nw=nw)
            constraint_thermal_limit_to_soft(pm, i, nw=nw)
        end

        for i in _PM.ids(pm, nw=nw, :busdc)                                               
            constraint_power_balance_dc_soft(pm, i, nw=nw)                                     
        end                                                                               
        for i in _PM.ids(pm, nw=nw, :branchdc)                                             
            constraint_ohms_dc_branch_soft(pm, i, nw=nw)                                
        end                                                                                
        for i in _PM.ids(pm, nw=nw, :convdc)                                               
            _PMACDC.constraint_converter_losses(pm, i, nw=nw)                                                          
            _PMACDC.constraint_conv_transformer(pm, i, nw=nw)                              
            _PMACDC.constraint_conv_reactor(pm, i, nw=nw)                                  
            _PMACDC.constraint_conv_filter(pm, i, nw=nw)
            constraint_converter_current(pm, i, nw=nw)                                    
            if pm.ref[:it][:pm][:nw][nw][:convdc][i]["islcc"] == 1            
                _PMACDC.constraint_conv_firing_angle(pm, i, nw=nw)                                    
            end
            constraint_pvdc_droop_control_smooth(pm, i, nw=nw)                                                                            
        end
 
    end
    objective_min_fuel_cost_scopf_soft(pm)
       
end



"""
This formulation is best used in conjunction with the contingency filters that find
violated contingencies in integrated ACDC grid.

"""

function run_scopf_soft_minlp(data, model_constructor, solver; kwargs...)
    # _PMACDC.process_additional_data!(data)
    return _PM.solve_model(data, model_constructor, solver, build_scopf_soft_minlp; ref_extensions = [_PM.ref_add_on_off_va_bounds!, _PMACDC.add_ref_dcgrid!, _PMSC.ref_c1!], multinetwork=true, kwargs...) 
end

function build_scopf_soft_minlp(pm::_PM.AbstractPowerModel)
    
    _PM.variable_bus_voltage(pm, nw=0)
    _PM.variable_gen_power(pm, nw=0)
    _PM.variable_branch_power(pm, nw=0)
    _PM.variable_branch_transform(pm, nw=0)
    _PM.constraint_model_voltage(pm, nw=0)

    _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=0)

    _PMACDC.variable_active_dcbranch_flow(pm, nw=0)       
    _PMACDC.variable_dcbranch_current(pm, nw=0)
    _PMACDC.variable_dc_converter(pm, nw=0)                        
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=0) 
    _PMACDC.constraint_voltage_dc(pm, nw=0)    
   
    variables_slacks(pm, nw=0)
    variable_converter_droop_binary(pm, nw=0)

    for i in _PM.ids(pm, nw=0, :ref_buses)                  
        _PM.constraint_theta_ref(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :bus)                        
        constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :branch)                     
        constraint_ohms_y_oltc_pst_from(pm, i, nw=0)
        constraint_ohms_y_oltc_pst_to(pm, i, nw=0)
        _PM.constraint_voltage_angle_difference(pm, i, nw=0)
        constraint_thermal_limit_from_soft(pm, i, nw=0)
        constraint_thermal_limit_to_soft(pm, i, nw=0)
    end

    for i in _PM.ids(pm, nw=0, :busdc)                                        
        constraint_power_balance_dc_soft(pm, i, nw=0)                                      
    end

    for i in _PM.ids(pm, nw=0, :branchdc)                                             
        constraint_ohms_dc_branch_soft(pm, i, nw=0)                                 
    end 

    for i in _PM.ids(pm, nw=0, :convdc)                                                
        _PMACDC.constraint_converter_losses(pm, i, nw=0)                                                            
        _PMACDC.constraint_conv_transformer(pm, i, nw=0)                               
        _PMACDC.constraint_conv_reactor(pm, i, nw=0)                                 
        _PMACDC.constraint_conv_filter(pm, i, nw=0)
        constraint_converter_current(pm, i, nw=0)                                     
        if pm.ref[:it][:pm][:nw][_PM.nw_id_default][:convdc][i]["islcc"] == 1                
            _PMACDC.constraint_conv_firing_angle(pm, i, nw=0)                                     
        end
        constraint_pvdc_droop_control_milp(pm, i, nw=0)                                                                        
    end                                                 

    contigency_ids = [id for id in _PM.nw_ids(pm) if id != 0]         
    
    for nw in contigency_ids

        _PM.variable_bus_voltage(pm, nw=nw, bounded=false)
        _PM.variable_gen_power(pm, nw=nw, bounded=false)
        _PM.variable_branch_power(pm, nw=nw)
        _PM.variable_branch_transform(pm, nw=nw)
        _PM.constraint_model_voltage(pm, nw=nw)

        _PMSC.variable_c1_shunt_admittance_imaginary(pm, nw=nw)
        _PMSC.variable_c1_response_delta(pm, nw=nw)

        _PMACDC.variable_active_dcbranch_flow(pm, nw=nw, bounded=false)       
        _PMACDC.variable_dcbranch_current(pm, nw=nw)
        _PMACDC.variable_dc_converter(pm, nw=nw)                        
        _PMACDC.variable_dcgrid_voltage_magnitude(pm, nw=nw)
        _PMACDC.constraint_voltage_dc(pm, nw=nw)
          
        variables_slacks(pm, nw=nw)
        variable_gen_response_binary(pm, nw=nw)
        variable_converter_droop_binary(pm, nw=nw)

        for i in _PM.ids(pm, nw=nw, :ref_buses)           
            _PM.constraint_theta_ref(pm, i, nw=nw)
        end

        gen_buses = _PM.ref(pm, nw=nw, :gen_buses)
        for i in _PM.ids(pm, nw=nw, :bus)                 
            constraint_power_balance_ac_shunt_dispatch_soft(pm, i, nw=nw)       

            # for a bus with active generators, fix voltage magnitude to base case until limits are reached, then switch PV/PQ bus
            if i in gen_buses
                constraint_gen_power_reactive_response_milp(pm, i, nw_1=0, nw_2=nw)
            end
        end

        response_gens = _PM.ref(pm, :response_gens, nw=nw)
        for (i,gen) in _PM.ref(pm, :gen, nw=nw)
            pg_base = _PM.var(pm, :pg, i, nw=0)

            # setup the linear response function until limits are reached or fix value to base case
            if i in response_gens
                constraint_gen_power_real_response_milp(pm, i, nw_1=0, nw_2=nw)
            else
                _PMSC.constraint_c1_gen_power_real_link(pm, i, nw_1=0, nw_2=nw)
            end
        end
  
        for i in _PM.ids(pm, nw=nw, :branch)                     
            constraint_ohms_y_oltc_pst_from(pm, i, nw=nw)
            constraint_ohms_y_oltc_pst_to(pm, i, nw=nw)
            _PM.constraint_voltage_angle_difference(pm, i, nw=nw)
            constraint_thermal_limit_from_soft(pm, i, nw=nw)
            constraint_thermal_limit_to_soft(pm, i, nw=nw)
        end

        for i in _PM.ids(pm, nw=nw, :busdc)                                               
            constraint_power_balance_dc_soft(pm, i, nw=nw)                                     
        end                                                                               
        for i in _PM.ids(pm, nw=nw, :branchdc)                                             
            constraint_ohms_dc_branch_soft(pm, i, nw=nw)                                
        end                                                                                
        for i in _PM.ids(pm, nw=nw, :convdc)                                               
            _PMACDC.constraint_converter_losses(pm, i, nw=nw)                                                          
            _PMACDC.constraint_conv_transformer(pm, i, nw=nw)                              
            _PMACDC.constraint_conv_reactor(pm, i, nw=nw)                                  
            _PMACDC.constraint_conv_filter(pm, i, nw=nw)
            constraint_converter_current(pm, i, nw=nw)                                    
            if pm.ref[:it][:pm][:nw][nw][:convdc][i]["islcc"] == 1            
                _PMACDC.constraint_conv_firing_angle(pm, i, nw=nw)                                    
            end
            constraint_pvdc_droop_control_milp(pm, i, nw=nw)                                                                            
        end
    end

    objective_min_fuel_cost_scopf_soft(pm)
end
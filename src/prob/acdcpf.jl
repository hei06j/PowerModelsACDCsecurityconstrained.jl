# ...existing code...
export solve_acdcpf

"""
    solve_acdcpf(file::String, model_type::Type, solver; kwargs...)

Parse an input file and solve an AC/DC power flow (ACDCPF) problem using the
specified PowerModels formulation and solver.

# Arguments
- `file::String` : Path to the input data file (e.g., MATPOWER `.m` file).
- `model_type::Type` : PowerModels model type (e.g., `ACPPowerModel`, `SOCBFPowerModel`).
- `solver` : JuMP solver object or solver factory (e.g., Ipopt).
- `kwargs...` : Optional keyword arguments forwarded to the underlying solver
  entry point (settings, reference extensions, etc.).

# Returns
- A solution dictionary (as returned by PowerModels) containing results such as
  objective value, variable values and solver termination status.

# Behavior
This function parses the provided file into a data dictionary, performs any
additional data processing (via `process_additional_data!`) and delegates to
`solve_acdcpf(data::Dict, ...)`. By default, reference extensions for DC grid,
PST, SSSC, flexible loads and DC generators are applied.
"""
function solve_acdcpf(file::String, model_type::Type, solver; kwargs...)
    data = _PM.parse_file(file)
    _PMACDC.process_additional_data!(data)
    return solve_acdcpf(data, model_type, solver; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.ref_add_pst!, _PMACDC.ref_add_sssc!, _PMACDC.ref_add_flex_load!, _PMACDC.ref_add_gendc!], kwargs...)
end

"""
    solve_acdcpf(data::Dict{String,Any}, model_type::Type, solver; kwargs...)

Solve an AC/DC power flow problem from an already-parsed data dictionary.

# Arguments
- `data::Dict{String,Any}` : Parsed network data dictionary (PowerModels format).
- `model_type::Type` : PowerModels model type to build the JuMP model.
- `solver` : JuMP solver object or factory.
- `kwargs...` : Forwarded keyword arguments (settings, ref_extensions override, etc.).

# Returns
- A solution dictionary produced by the PowerModels solve pipeline.

# Behavior
Builds and solves the problem using the power model solve entry point (`_PM.solve_model`).
This wrapper applies the same set of default reference extensions as the file-based
entrypoint. Use `ref_extensions` in `kwargs` to override or add additional references.
"""
function solve_acdcpf(data::Dict{String,Any}, model_type::Type, solver; kwargs...)
    return _PM.solve_model(data, model_type, solver, build_acdcpf; ref_extensions = [_PMACDC.add_ref_dcgrid!, _PMACDC.ref_add_pst!, _PMACDC.ref_add_sssc!, _PMACDC.ref_add_flex_load!, _PMACDC.ref_add_gendc!], kwargs...)
end

"""
    build_acdcpf(pm::_PM.AbstractPowerModel)

Builds the JuMP model for an AC/DC power flow (ACDCPF) formulation.

# Arguments
- `pm::_PM.AbstractPowerModel` : PowerModels internal model holder (pm).

# Details
Creates variables and constraints appropriate for non-optimization power flow
(bounded = false variants are used). The builder:
- Adds bus voltage, generator active power, branch power and storage power variables.
- For SOCBF formulations, also adds branch current variables.
- Adds DC-specific variables (active DC branch flows, DC branch current, converters,
  DC grid voltage magnitudes, DC generators).
- Adds flexible and fixed load handling, PST and SSSC variables when present.
- Adds model constraints including:
  - voltage model constraints
  - DC voltage constraints
  - reference bus constraints (theta reference and voltage setpoints)
  - AC power balance and PV bus handling (generator setpoints / voltage setpoints)
  - branch equations (Ohm's law or SOCBF-specific constraints)
  - DC power balance and DC branch Ohm's law
  - converter-specific constraints (transformer, reactor, filter, setpoints,
    losses and current limits)
- For converters with droop or DC voltage setpoint control, appropriate constraints
  are created depending on converter `type_dc` and `type_ac`.
"""
function build_acdcpf(pm::_PM.AbstractPowerModel)
    _PM.variable_bus_voltage(pm, bounded = false)
    _PM.variable_gen_power(pm, bounded = false)
    _PM.variable_branch_power(pm, bounded = false)
    _PM.variable_storage_power(pm, bounded = false)

    # dirty, should be improved in the future TODO
    if typeof(pm) <: _PM.SOCBFPowerModel
        _PM.variable_branch_current(pm, bounded = false)
    end

    _PMACDC.variable_active_dcbranch_flow(pm, bounded = false)
    _PMACDC.variable_dcbranch_current(pm, bounded = false)
    _PMACDC.variable_dc_converter(pm, bounded = false)
    _PMACDC.variable_dcgrid_voltage_magnitude(pm, bounded = false)
    _PMACDC.variable_dcgenerator_power(pm; bounded = false)
    _PMACDC.variable_flexible_demand(pm, bounded = false)
    _PMACDC.variable_pst(pm, bounded = false)
    _PMACDC.variable_sssc(pm, bounded = false)

    _PM.constraint_model_voltage(pm)
    _PMACDC.constraint_voltage_dc(pm)


    for (i,bus) in _PM.ref(pm, :ref_buses)
        @assert bus["bus_type"] == 3
        _PM.constraint_theta_ref(pm, i)
        _PM.constraint_voltage_magnitude_setpoint(pm, i)
    end

    for (i, bus) in _PM.ref(pm, :bus)# _PM.ids(pm, :bus)
        _PMACDC.constraint_power_balance_ac(pm, i)
        # PV Bus Constraints
        if length(_PM.ref(pm, :bus_gens, i)) > 0 && !(i in _PM.ids(pm,:ref_buses))
            for j in _PM.ref(pm, :bus_gens, i)
                _PM.constraint_gen_setpoint_active(pm, j)
                if  bus["bus_type"] == 2
                    _PM.constraint_voltage_magnitude_setpoint(pm, i)
                elseif bus["bus_type"] == 1
                    _PM.constraint_gen_setpoint_active(pm, j)
                end
            end
        end
    end

    for i in _PM.ids(pm, :branch)
        # dirty, should be improved in the future TODO
        if typeof(pm) <: _PM.SOCBFPowerModel
            _PM.constraint_power_losses(pm, i)
            _PM.constraint_voltage_magnitude_difference(pm, i)
            _PM.constraint_branch_current(pm, i)
        else
            _PM.constraint_ohms_yt_from(pm, i)
            _PM.constraint_ohms_yt_to(pm, i)
        end
    end
    for i in _PM.ids(pm, :flex_load)
        _PMACDC.constraint_total_flexible_demand(pm, i)
    end
    
    for i in _PM.ids(pm, :fixed_load) 
        _PMACDC.constraint_total_fixed_demand(pm, i)
    end

    for i in _PM.ids(pm, :busdc)
        _PMACDC.constraint_power_balance_dc(pm, i)
    end
    for i in _PM.ids(pm, :branchdc)
        _PMACDC.constraint_ohms_dc_branch(pm, i)
    end

    if !isempty(_PM.ids(pm, :gendc)) 
        for i in _PM.ids(pm, :gendc)
            _PMACDC.constraint_dcgenerator_voltage_and_power(pm, i)
        end
    end

    for (c, conv) in _PM.ref(pm, :convdc)
        # _PMACDC.constraint_conv_transformer(pm, c)
        # _PMACDC.constraint_conv_reactor(pm, c)
        _PMACDC.constraint_conv_filter(pm, c)
        constraint_conv_transformer(pm, c)
        constraint_conv_reactor(pm, c)
        if conv["type_dc"] == 2
            _PMACDC.constraint_dc_voltage_magnitude_setpoint(pm, c)
        elseif conv["type_dc"] == 3 || conv["type_dc"] == 4
            if typeof(pm) <: _PM.AbstractACPModel || typeof(pm) <: _PM.AbstractACRModel
                _PMACDC.constraint_dc_droop_control(pm, c)
            else
                Memento.warn(_PM._LOGGER, join(["Droop only defined for ACP and ACR formulations, converter ", c, " will be treated as type 2"]))
                _PMACDC.constraint_dc_voltage_magnitude_setpoint(pm, c)
            end
        else
            _PMACDC.constraint_active_conv_setpoint(pm, c)
        end
        if conv["type_ac"] == 2
            if haskey(conv, "acq_droop") && conv["acq_droop"] == 1 # AC voltage droop control
                _PMACDC.constraint_ac_voltage_droop_control(pm, c)
            else # Constant AC voltage control
                _PMACDC.constraint_voltage_magnitude_setpoint(pm, conv["busac_i"])
            end
        else
            _PMACDC.constraint_reactive_conv_setpoint(pm, c)
        end
        _PMACDC.constraint_converter_losses(pm, c)
        # _PMACDC.constraint_converter_current(pm, c)
        constraint_converter_current(pm, c)
    end
end
# ...existing code...
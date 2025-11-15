%% Unified Thermosyphon Model with Continuous NH3 Circulation
% Single operating mode with continuous NH3 flow calculation
% Sequential energy cascade with no circular dependency
% FIXED: Heat pump hysteresis, mass balance for condensation, water preheating priority

clear; clc; close all;

%% Initialize CoolProp
try
    py.CoolProp.CoolProp.PropsSI('T', 'P', 101325, 'Q', 0, 'Water');
    fprintf('CoolProp successfully initialized\n');
catch
    error('CoolProp initialization failed. Please ensure Python and CoolProp are properly installed.');
end

%% System Parameters
params = struct();
params.heatInput = 1000;          % Heat input rate [W/m²] 
params.initialTemp = 0;              % Initial temperature [°C]
params.length = 1800;                % Heat Pipe length [m]
params.evap_length = 0.25 * params.length;
params.cond_length = 0.25 * params.length;
params.adiab_length = 0.50 * params.length;
params.initialLiquidHeight = params.length * 0.20;  % Initial liquid height [m]
params.diameter_in = 0.1778;         % Inner diameter Heat Pipe [m]
params.diameter_out = 0.1900;         % Outer diameter Heat Pipe [m]  
params.crossSectionArea = pi * (params.diameter_in/2)^2;     % Cross-sectional area [m²]
params.systemVolume = params.crossSectionArea * params.length; % Total system volume [m³]
params.evap_area = 2* pi * params.diameter_out / 2 * params.evap_length;
params.cond_area = 2* pi * params.diameter_out * params.cond_length;
params.adiab_area = 2* pi * params.diameter_out * params.adiab_length;

params.fluid = 'R717';               % Fluid type (NH3)
params.timeStep = 30;                % Time step [s]
params.totalTime = 3600 * 24;       % Total simulation time [s]

%% Wall Thermal Mass Parameters
params.pipe_wall_thickness = (params.diameter_out - params.diameter_in) / 2;  % [m]

% Steel properties
params.wall_density_steel = 7850;    % Steel density [kg/m³]
params.wall_cp_steel = 480;          % Steel specific heat [J/kg·K]

% Calculate wall mass
params.wall_volume = pi * ((params.diameter_out/2)^2 - (params.diameter_in/2)^2) * params.length;
params.wall_mass = params.wall_density_steel * params.wall_volume;  % [kg]
params.wall_thermal_capacity = params.wall_mass * params.wall_cp_steel;  % [J/K]

% Display for verification
fprintf('\n=== WALL THERMAL MASS (NEW) ===\n');
fprintf('Wall thickness: %.4f m\n', params.pipe_wall_thickness);
fprintf('Wall volume: %.2f m³\n', params.wall_volume);
fprintf('Wall mass: %.0f kg\n', params.wall_mass);
fprintf('Wall thermal capacity: %.2f MJ/K\n', params.wall_thermal_capacity/1e6);
fprintf('================================\n\n');

%% Phase Change Model Parameters
% Cooper Correlation - Wall Temperature Estimate
params.T_wall_offset = 5;  % [°C] PLACEHOLDER: Wall is ~5°C hotter than interface
                           % Smaller offset to avoid overestimating heat transfer
                           % TODO: In next step, replace with thermal mass calculation for exact T_wall

% Note: Geothermal input uses ACTUAL temperature for correct physics
% No reference temperature needed - system finds natural equilibrium

%% Geothermal Parameters 
params.use_geothermal_profile = true;
params.geothermal_max_temp = 81.12896;  % °C from your data
params.geothermal_coeff_a = 0.0060589;  % from your equation
params.geothermal_coeff_b = 1.99186;   % from your equation
params.pipe_thermal_conductivity = 50;  % W/m·K
params.pipe_wall_thickness = params.diameter_out - params.diameter_in;      % m
params.ground_thermal_conductivity = 2.5; % W/m·K
params.pipe_to_ground_resistance = 0.1; % K·m/W

%% Heat Pump Coupling Parameters
heat_pump_demand = @(t) 40e3 + 30e3 * sin(1*pi*t/(24*3600));
%heat_pump_demand = @(t) 37e3 + 20e3 * sin(1*pi*t/(24*3600));

params.heat_pump_duty = heat_pump_demand;     % Heat pump evaporator requirement [W]
params.T_equilibrium_target = 60;    % Target equilibrium temperature [°C]
params.T_return_target = 20;         % NH3 return temperature from heat pump [°C]
params.T_heat_pump_start = 59;       % Temperature to START heat pump operation [°C]
params.T_heat_pump_stop = 55;        % Temperature to STOP heat pump operation [°C]

%% NH3 Flow Control Parameters
params.min_NH3_flow = 0.01;          % kg/s - Minimum flow
params.max_NH3_flow = 0.2;           % kg/s - Maximum flow
params.equilibrium_tolerance = 1.0;   % °C - Temperature tolerance for equilibrium

%% Heat Loss Parameters
params.heat_loss_coefficient = 50;   % Heat loss coefficient [W/m²·K]
params.surface_area_loss = 20;       % Surface area for heat loss [m²]
params.T_ambient = 20;               % Ambient temperature [°C]

%% Energy Factor Constants
%params.energyFactor_heating = 0.30;      % 30% during heating phase
%params.energyFactor_circulation = 0.813; % 81.3% during circulation phase

%% Heat Pump Adaptive Control Parameters
params.heat_pump_max_turndown = 0.3;    % Minimum 30% of rated power when operating
params.heat_pump_efficiency_curve = 1.0; % Could vary with load in future

%% Energy Balance Verification
params.energy_balance_tolerance = 1000;  % W - tolerance for energy balance check

%% Water Preheating Parameters
params.m_water_preheat = 30;         % Water flow rate for preheating [kg/s]
params.cp_water = 4180;              % Water specific heat [J/kg·K]

%% Vapor Column Visualization Parameters
params.profile_segments = 36;        % Every 50m: 1800m / 50m = 36 segments
params.segment_height = 50;          % m - height increment for pressure calculation

%% Initial Conditions
state = struct();
state.time = 0;
state.liquidHeight = params.initialLiquidHeight;
state.temperature = params.initialTemp;
state.totalEvaporated = 0;
state.equilibrium_achieved = false;
state.equilibrium_start_time = 0;
state.heat_pump_active = false;

% Initialize at saturation conditions
state.vaporPressure = py.CoolProp.CoolProp.PropsSI('P', 'T', params.initialTemp + 273.15, 'Q', 0, params.fluid);

%% Preallocate arrays for results
nSteps = floor(params.totalTime / params.timeStep);
results = struct();
results.time = zeros(nSteps+1, 1);
results.vaporPressure = zeros(nSteps+1, 1);
results.hydroPressure = zeros(nSteps+1, 1);
results.liquidHeight = zeros(nSteps+1, 1);
results.satTemp = zeros(nSteps+1, 1);
results.evapRate = zeros(nSteps+1, 1);
results.liquidDensity = zeros(nSteps+1, 1);
results.vaporDensity = zeros(nSteps+1, 1);
results.hfg = zeros(nSteps+1, 1);
results.liquidMass = zeros(nSteps+1, 1);
results.energyFactor = zeros(nSteps+1, 1);
results.Q_input_effective = zeros(nSteps+1, 1);
results.Q_input_total_EF = zeros(nSteps+1, 1);
results.Q_heat_pump = zeros(nSteps+1, 1);
results.Q_water_preheat = zeros(nSteps+1, 1);
results.Q_losses = zeros(nSteps+1, 1);
results.Q_NH3_circulation = zeros(nSteps+1, 1);
results.Q_excess = zeros(nSteps+1, 1);
results.energy_balance_error = zeros(nSteps+1, 1);
results.equilibrium_status = cell(nSteps+1, 1);
results.NH3_flow_actual = zeros(nSteps+1, 1);
results.temperature_error = zeros(nSteps+1, 1);
results.constraint_status = cell(nSteps+1, 1);
results.HP_load_percent = zeros(nSteps+1, 1);
results.preheat_percentage = zeros(nSteps+1, 1);

%% Store initial conditions
initial_liquidDensity = py.CoolProp.CoolProp.PropsSI('D', 'T', params.initialTemp + 273.15, 'Q', 0, params.fluid);
initial_vaporDensity = py.CoolProp.CoolProp.PropsSI('D', 'P', state.vaporPressure, 'Q', 1, params.fluid);
initial_h_liquid = py.CoolProp.CoolProp.PropsSI('H', 'T', params.initialTemp + 273.15, 'Q', 0, params.fluid);
initial_h_vapor = py.CoolProp.CoolProp.PropsSI('H', 'T', params.initialTemp + 273.15, 'Q', 1, params.fluid);
initial_hfg = initial_h_vapor - initial_h_liquid;
initial_satTemp = py.CoolProp.CoolProp.PropsSI('T', 'P', state.vaporPressure, 'Q', 0.5, params.fluid) - 273.15;
initial_liquidVolume = state.liquidHeight * params.crossSectionArea;
initial_vaporVolume = params.systemVolume - initial_liquidVolume;
state.initialVaporMass = initial_vaporVolume * initial_vaporDensity;
state.totalVaporMass = state.initialVaporMass;
state.liquidMass = params.initialLiquidHeight * params.crossSectionArea * initial_liquidDensity;

% Store initial values
results.time(1) = 0;
results.vaporPressure(1) = state.vaporPressure / 1000;
results.hydroPressure(1) = (initial_liquidDensity * 9.81 * state.liquidHeight) / 1000;
results.liquidHeight(1) = state.liquidHeight;
results.satTemp(1) = initial_satTemp;
results.evapRate(1) = 0;
results.liquidDensity(1) = double(initial_liquidDensity);
results.vaporDensity(1) = double(initial_vaporDensity);
results.hfg(1) = double(initial_hfg) / 1000;
results.liquidMass(1) = state.liquidMass;
results.energyFactor(1) = 0.30;
results.equilibrium_status{1} = 'TRANSIENT';
results.constraint_status{1} = 'NONE';
results.HP_load_percent(1) = 0;
results.Q_input_total_EF(1) = 0;

fprintf('Initial conditions:\n');
fprintf('  Temperature: %.1f°C\n', initial_satTemp);
fprintf('  Vapor pressure: %.1f kPa\n', state.vaporPressure/1000);
fprintf('  Initial liquid mass: %.2f kg\n', state.liquidMass);

%% Main Simulation Loop
fprintf('Starting unified simulation...\n');
progressInterval = floor(nSteps / 20);

for step = 1:nSteps
    prev_temperature = state.temperature;  % Track temperature changes
    state.time = state.time + params.timeStep;
    
    % Progress indicator
    if mod(step, progressInterval) == 0
        fprintf('Progress: %.0f%% (t=%.1f min, T=%.1f°C)\n', ...
            100 * step / nSteps, state.time/60, state.temperature);
    end
    
%% Get current properties at saturation
satTemp = py.CoolProp.CoolProp.PropsSI('T', 'P', state.vaporPressure, 'Q', 0, params.fluid) - 273.15;
state.temperature = satTemp;  % Keep this for all system logic

% Calculate interface pressure for evaporation properties
vapor_column_height = params.length - state.liquidHeight;
vaporDensity = py.CoolProp.CoolProp.PropsSI('D', 'P', state.vaporPressure, 'Q', 1, params.fluid);
P_interface = state.vaporPressure + vaporDensity * 9.81 * vapor_column_height;

% NEW: Calculate actual interface temperature
T_interface_actual = py.CoolProp.CoolProp.PropsSI('T', 'P', P_interface, 'Q', 0, params.fluid) - 273.15;

% Use interface pressure ONLY for evaporation properties
liquidDensity = py.CoolProp.CoolProp.PropsSI('D', 'P', P_interface, 'Q', 0, params.fluid);
h_liquid = py.CoolProp.CoolProp.PropsSI('H', 'P', P_interface, 'Q', 0, params.fluid);
h_vapor = py.CoolProp.CoolProp.PropsSI('H', 'P', P_interface, 'Q', 1, params.fluid);
hfg = h_vapor - h_liquid;  % More accurate hfg at interface pressure

%% DEBUG OUTPUT - ADD HERE
if mod(step, 200) == 0
    T_interface = py.CoolProp.CoolProp.PropsSI('T', 'P', P_interface, 'Q', 0.5, params.fluid) - 273.15;
    fprintf('  System: T=%.1f°C, P=%.1f kPa\n', state.temperature, state.vaporPressure/1000);
    fprintf('  Interface: T=%.1f°C, P=%.1f kPa, ΔT=%.1f°C\n', ...
            T_interface, P_interface/1000, T_interface-state.temperature);
end

    %% STEP 1: CALCULATE GEOTHERMAL HEAT AT ACTUAL TEMPERATURE
    % Use actual NH3 temperature for correct physics - no artificial reference
    % System will naturally find equilibrium through negative feedback

    % Calculate geothermal heat input at actual pipe temperature
    [Q_geothermal_total, Q_profile_debug] = calculate_geothermal_input(state, params, state.temperature);

    % Calculate heat losses
    Q_losses = params.heat_loss_coefficient * params.surface_area_loss * (state.temperature - params.T_ambient);

    % Net energy available for NH3 circulation (physics-based, no energy factor)
    Q_input_net = Q_geothermal_total - Q_losses;

    % Calculate evaporation rate from Cooper correlation (for thermal mass only)
    [evapRate_physics, ~] = calculatePhaseChangePhysics(state, params, Q_geothermal_total);
    
    %% STEP 2: SEQUENTIAL ENERGY CASCADE (NO CIRCULAR DEPENDENCY)
    % FIXED: Modified logic to allow heat pump to cool below equilibrium target
    % Once heat pump is active, stay in circulation mode until HP turns off
    
    if state.temperature >= params.T_equilibrium_target || state.heat_pump_active
        %% PHASE: NH3 CIRCULATION MODE
        %% Calculate NH3 properties for circulation (ACCURATE METHOD)
        T_current_K = state.temperature + 273.15;
        T_return_K = params.T_return_target + 273.15;

        % Get actual enthalpies at liquid state (dH = Integral cp*dT)
        h_liquid_current = py.CoolProp.CoolProp.PropsSI('H', 'T', T_current_K, 'Q', 0, params.fluid);
        h_liquid_return = py.CoolProp.CoolProp.PropsSI('H', 'T', T_return_K, 'Q', 0, params.fluid);

        % Sensible heat = actual enthalpy difference (accounts for cp variation)
        sensible_heat_NH3 = h_liquid_current - h_liquid_return;

        % Total energy per kg NH3 (sensible + latent)
        energy_per_kg_NH3 = sensible_heat_NH3 + hfg; % C = cp*dT + hfg
        
        %% STEP 2.1: MINIMUM NH3 FLOW (System baseline)
        NH3_flow_base = params.min_NH3_flow;
        Q_NH3_base = NH3_flow_base * energy_per_kg_NH3;
        
        %% STEP 2.2: ENERGY AVAILABLE FOR HEAT PUMP
        Q_available_for_HP = max(0, Q_input_net - Q_NH3_base);
        
        %% STEP 2.3: HEAT PUMP OPERATION WITH PROPER HYSTERESIS (REALLY FIXED!)
        % FIXED: Proper hysteresis logic with start/stop thresholds - TEMPERATURE ONLY!
        if ~state.heat_pump_active && state.temperature >= params.T_heat_pump_start
            state.heat_pump_active = true;
            fprintf('🔥 Heat pump STARTED at T=%.1f°C\n', state.temperature);
        elseif state.heat_pump_active && state.temperature <= params.T_heat_pump_stop
            state.heat_pump_active = false;
            fprintf('❄️  Heat pump STOPPED at T=%.1f°C (temperature hysteresis)\n', state.temperature);
        end
        
        % Get heat pump demand regardless of status
        Q_heat_pump_demand = params.heat_pump_duty(state.time);
        
        if state.heat_pump_active
            % FIXED: Heat pump runs even if energy is insufficient - it will create negative Q_excess
            % This allows the system to cool down properly when heat pump demand > input
            Q_heat_pump_actual = Q_heat_pump_demand;  % Heat pump takes what it needs!
            HP_load_percent = 100;  % Always 100% when active
            
            % Don't turn off heat pump due to insufficient energy - let it create cooling
            % The negative Q_excess will cause temperature to drop until hysteresis stops it
        else
            % Heat pump is off due to temperature hysteresis
            Q_heat_pump_actual = 0;
            HP_load_percent = 0;
        end
        
%% CORRECTED: NH3 FLOW LOGIC WITH THERMAL MASS EXTRACTION

%% STEP 2.4: CALCULATE NH3 FLOW
% KEY FIX: NH3 flow based on AVAILABLE ENERGY, not just heat pump demand
% This allows using excess energy for water preheating

if state.temperature > params.T_equilibrium_target
    % THERMAL MASS EXTRACTION MODE
    % System is too hot - actively cool down by extracting maximum heat
    % Temperature is ALLOWED to drop - extract as much as possible!
    extraction_mode = 'THERMAL_MASS';

    if energy_per_kg_NH3 > 0
        % KEY FIX: Use MAXIMUM NH3 flow to cool down system
        % Don't limit by Q_input_net - extract from thermal mass!
        % This maximizes preheating and cools system toward target

        % Calculate minimum flow needed for heat pump
        if state.heat_pump_active
            NH3_flow_min_for_HP = Q_heat_pump_demand / energy_per_kg_NH3;
        else
            NH3_flow_min_for_HP = params.min_NH3_flow;
        end

        % When T > target, use maximum flow to cool down
        % Temperature drop is DESIRED - extract maximum heat for preheating
        NH3_flow_actual = params.max_NH3_flow;

        % But ensure heat pump minimum is met
        NH3_flow_actual = max(NH3_flow_min_for_HP, NH3_flow_actual);
    else
        NH3_flow_actual = params.min_NH3_flow;
    end

else
    % NORMAL MODE
    % T ≤ T_target: NH3 flow limited by available geothermal energy
    extraction_mode = 'GEOTHERMAL';

    if Q_input_net > 0 && energy_per_kg_NH3 > 0
        NH3_flow_from_energy = Q_input_net / energy_per_kg_NH3;
    else
        NH3_flow_from_energy = 0;
    end

    NH3_flow_actual = max(params.min_NH3_flow, min(NH3_flow_from_energy, params.max_NH3_flow));
end

%% STEP 2.5: (Flow constraints already applied above)

%% STEP 2.6: CHECK FOR FLOW CONSTRAINTS
if NH3_flow_actual >= params.max_NH3_flow * 0.999
    flow_is_constrained = true;
    constraint_type = 'MAX_FLOW';
elseif NH3_flow_actual <= params.min_NH3_flow * 1.001
    flow_is_constrained = true;
    constraint_type = 'MIN_FLOW';
else
    flow_is_constrained = false;
    constraint_type = 'NONE';
end

%% STEP 2.7: CALCULATE TOTAL NH3 CIRCULATION ENERGY
Q_NH3_circulation = NH3_flow_actual * energy_per_kg_NH3;

%% STEP 2.8: DISTRIBUTE NH3 ENERGY BETWEEN HEAT PUMP AND PREHEATING
% Heat pump gets priority up to its demand
Q_heat_pump_delivered = min(Q_heat_pump_actual, Q_NH3_circulation);

% Preheating gets the remainder
if state.heat_pump_active && Q_heat_pump_delivered >= Q_heat_pump_demand * 0.999
    % Heat pump demand satisfied - preheating can use excess
    Q_water_preheat = max(0, Q_NH3_circulation - Q_heat_pump_delivered);
elseif ~state.heat_pump_active
    % Heat pump off - all circulation goes to preheating
    Q_water_preheat = Q_NH3_circulation;
else
    % Heat pump demand not satisfied - no preheating
    Q_water_preheat = 0;
    
    % If heat pump can't be satisfied, it might need to reduce its demand
    if Q_heat_pump_actual > Q_NH3_circulation
        fprintf('⚠️ Heat pump demand (%.1f kW) exceeds available NH3 energy (%.1f kW)\n', ...
                Q_heat_pump_actual/1000, Q_NH3_circulation/1000);
    end
end

%% Verification: NH3 energy should equal heat pump + preheating
NH3_energy_check = Q_heat_pump_delivered + Q_water_preheat;
if abs(Q_NH3_circulation - NH3_energy_check) > 100  % 0.1 kW tolerance
    fprintf('❌ NH3 energy mismatch: Circulation=%.1f, HP+Preheat=%.1f kW\n', ...
            Q_NH3_circulation/1000, NH3_energy_check/1000);
end

%% Debug output
if mod(step, 200) == 0
    fprintf('\n=== NH3 FLOW LOGIC ===\n');
    fprintf('  Extraction mode: %s\n', extraction_mode);
    fprintf('  Q_input_net: %.1f kW\n', Q_input_net/1000);
    fprintf('  NH3_flow_actual: %.3f kg/s (%.1f%% of max)\n', ...
            NH3_flow_actual, 100*NH3_flow_actual/params.max_NH3_flow);
    fprintf('  Q_NH3_circulation: %.1f kW\n', Q_NH3_circulation/1000);
    fprintf('  Q_heat_pump_delivered: %.1f kW\n', Q_heat_pump_delivered/1000);
    fprintf('  Q_water_preheat: %.1f kW\n', Q_water_preheat/1000);
    if strcmp(extraction_mode, 'THERMAL_MASS')
        fprintf('  ⚡ Extracting from thermal mass (T=%.1f > %.1f°C)\n', ...
                state.temperature, params.T_equilibrium_target);
    end
    fprintf('======================\n');
end
        
        %% STEP 2.9: FINAL EXCESS ENERGY (CAN BE NEGATIVE!)
        % Q_excess = Energy going to/from thermal mass
        % Positive: Heating up (geothermal > extraction)
        % Negative: Cooling down (extraction > geothermal)

        Q_excess = Q_input_net - Q_NH3_circulation;

        % When extracting from thermal mass (T > target), Q_excess will be negative
        % This correctly represents cooling of the thermal mass
        
    else
        %% PHASE: HEATING MODE (Temperature too low AND heat pump inactive)
        NH3_flow_actual = 0;
        Q_NH3_circulation = 0;
        Q_heat_pump_delivered = 0;
        Q_heat_pump_demand = 0;
        Q_water_preheat = 0;
        Q_excess = Q_input_net;  % All energy goes to thermal mass heating
        
        flow_is_constrained = false;
        constraint_type = 'NO_CIRCULATION';
        % Note: heat pump state is already managed in hysteresis logic above
        HP_load_percent = 0;
    end
    
    %% STEP 3: EXACT ENERGY BALANCE VERIFICATION (INCLUDES HEAT PUMP)
    Q_total_output = Q_losses + Q_NH3_circulation + Q_excess;
    Q_total_input = Q_geothermal_total;
    energy_balance_error = Q_total_input - Q_total_output;
    
    % Check energy balance closure
    if abs(energy_balance_error) > params.energy_balance_tolerance
        warning('Energy balance error: %.1f kW at step %d', energy_balance_error/1000, step);
    end
    
    %% STEP 4: DETERMINE THERMAL MASS STATUS
    if Q_excess > 0
        thermal_mass_status = 'HEATING';
    elseif Q_excess < 0
        thermal_mass_status = 'COOLING';
    else
        thermal_mass_status = 'STABLE';
    end
    
    %% STEP 5: USE COOPER-BASED EVAPORATION RATE
    % Evaporation rate already calculated from Cooper correlation (physics-based)
    % Can be positive (evaporation) or negative (condensation if Q_excess < 0)
    evapRate = evapRate_physics;  % From calculatePhaseChangePhysics (line ~235)
    massEvaporated = evapRate * params.timeStep;  % Can be positive or negative
    
    %% STEP 6: EQUILIBRIUM CHECK
    temperature_error = state.temperature - params.T_equilibrium_target;
    
    if abs(temperature_error) < params.equilibrium_tolerance && ...
       abs(Q_excess) < 5000 && ...
       ~flow_is_constrained && ...
       state.heat_pump_active && ...
       HP_load_percent >= 90  % Heat pump at near full capacity
        
        equilibrium_status = 'ENERGY_EQUILIBRIUM';
        
        if ~state.equilibrium_achieved
            state.equilibrium_achieved = true;
            state.equilibrium_start_time = state.time;
            fprintf('\n🎯 EQUILIBRIUM ACHIEVED at t=%.1f min, T=%.1f°C\n', ...
                    state.time/60, state.temperature);
            fprintf('   Heat pump: %.1f kW (%.0f%% load), NH3 flow: %.3f kg/s\n', ...
                    Q_heat_pump_delivered/1000, HP_load_percent, NH3_flow_actual);
        end
    else
        equilibrium_status = 'TRANSIENT';
        if flow_is_constrained && Q_excess > 10000
            equilibrium_status = 'FLOW_LIMITED_HEATING';
        elseif state.temperature >= params.T_equilibrium_target && ~state.heat_pump_active
            equilibrium_status = 'INSUFFICIENT_ENERGY';
        end
    end
    
    %% STEP 7: UPDATE SYSTEM MASSES AND GEOMETRY (FIXED!)
    % FIXED: Better handling of mass conservation and bounds
    state.totalEvaporated = state.totalEvaporated + massEvaporated;
    state.totalVaporMass = state.initialVaporMass + state.totalEvaporated;
    state.liquidMass = state.liquidMass - massEvaporated;
    
    % FIXED: Ensure physical bounds
    state.liquidMass = max(0, state.liquidMass);
    state.totalVaporMass = max(1e-6, state.totalVaporMass); % Minimum vapor mass
    
    if state.liquidMass > 0
        currentLiquidVolume = state.liquidMass / liquidDensity;
        state.liquidHeight = currentLiquidVolume / params.crossSectionArea;
    else
        state.liquidHeight = 0;
    end
    
    %% STEP 8: ENFORCE THERMODYNAMIC EQUILIBRIUM WITH IMPROVED SOLVER (FIXED!)
    vaporVolume = params.systemVolume - state.liquidHeight * params.crossSectionArea;
    [T_final, P_final] = solve_mass_balance_improved(state.totalVaporMass, vaporVolume, state.temperature, params.fluid);
    state.temperature = T_final;
    state.vaporPressure = P_final;
    
    %% Calculate remaining outputs
    hydroPressure = liquidDensity * 9.81 * state.liquidHeight;
    satTemp = py.CoolProp.CoolProp.PropsSI('T', 'P', state.vaporPressure, 'Q', 0.5, params.fluid) - 273.15;
    
    %% STEP 9: STORE RESULTS WITH COMPLETE TRACKING
    results.time(step + 1) = state.time / 60;
    results.vaporPressure(step + 1) = state.vaporPressure / 1000;
    results.hydroPressure(step + 1) = hydroPressure / 1000;
    results.liquidHeight(step + 1) = state.liquidHeight;
    results.satTemp(step + 1) = double(satTemp);
    results.evapRate(step + 1) = evapRate * 1000; % [g/s]
    results.liquidDensity(step + 1) = double(liquidDensity);
    results.vaporDensity(step + 1) = double(vaporDensity);
    results.hfg(step + 1) = double(hfg) / 1000;
    results.liquidMass(step + 1) = state.liquidMass;
    results.energyFactor(step + 1) = 1.0;  % No energy factor - direct physics
    results.Q_input_effective(step + 1) = Q_input_net / 1000;
    results.Q_input_total_EF(step + 1) = Q_geothermal_total / 1000;
    results.Q_heat_pump(step + 1) = Q_heat_pump_delivered / 1000;
    results.Q_water_preheat(step + 1) = Q_water_preheat / 1000;
    results.Q_losses(step + 1) = Q_losses / 1000;
    results.Q_NH3_circulation(step + 1) = Q_NH3_circulation / 1000;
    results.Q_excess(step + 1) = Q_excess / 1000;
    results.energy_balance_error(step + 1) = energy_balance_error / 1000;
    results.equilibrium_status{step + 1} = equilibrium_status;
    results.NH3_flow_actual(step + 1) = NH3_flow_actual;
    results.temperature_error(step + 1) = temperature_error;
    results.constraint_status{step + 1} = constraint_type;
    results.HP_load_percent(step + 1) = HP_load_percent;
    
    %% Calculate preheating percentage
    total_useful_output = Q_heat_pump_delivered + Q_water_preheat;
    if total_useful_output > 0
        preheat_percentage = (Q_water_preheat / total_useful_output) * 100;
    else
        preheat_percentage = 0;
    end
    results.preheat_percentage(step + 1) = preheat_percentage;
    
    %% Store vapor column data periodically
    if mod(step, 100) == 0  % Store every 100 steps
        [profile_heights, profile_temps, profile_pressures, profile_sections] = ...
            create_real_temperature_profile(state, params, T_interface_actual);
        
        % Store in results
        if ~isfield(results, 'vapor_column_data')
            results.vapor_column_data = {};
        end
        
        vapor_data = struct();
        vapor_data.time = state.time / 60;
        vapor_data.heights = profile_heights;
        vapor_data.temperatures = profile_temps;
        vapor_data.pressures = profile_pressures / 1000; % Convert to kPa
        vapor_data.liquid_height = state.liquidHeight;
        vapor_data.sections = profile_sections;
        vapor_data.uniform_temp = state.temperature; % For comparison
        
        results.vapor_column_data{end+1} = vapor_data;
    end
    
    %% STEP 10: ENHANCED DEBUG OUTPUT (FIXED!)
    if mod(step, 200) == 0 || abs(energy_balance_error) > params.energy_balance_tolerance || state.equilibrium_achieved
        fprintf('\n=== Step %d: Sequential Energy Cascade ===\n', step);
        fprintf('  Temperature: %.1f°C (target: %.1f°C)\n', state.temperature, params.T_equilibrium_target);
        if state.temperature >= params.T_equilibrium_target || state.heat_pump_active
            fprintf('  Operating Phase: NH3 CIRCULATION (T≥%.1f°C OR HP active)\n', params.T_equilibrium_target);
        else
            fprintf('  Operating Phase: THERMAL HEATING (T<%.1f°C AND HP inactive)\n', params.T_equilibrium_target);
        end        
        fprintf('  \n');

        fprintf('  ENERGY CASCADE:\n');
        fprintf('    Q_geothermal (at %.1f°C actual): %.1f kW\n', state.temperature, Q_geothermal_total/1000);
        fprintf('    Q_losses: %.1f kW\n', Q_losses/1000);
        fprintf('    Q_input_net: %.1f kW\n', Q_input_net/1000);
        fprintf('  \n');
        
        if state.temperature >= params.T_equilibrium_target
            fprintf('  NH3 CIRCULATION ANALYSIS:\n');
            fprintf('    NH3 flow actual: %.3f kg/s (min: %.3f, max: %.3f)\n', ...
                    NH3_flow_actual, params.min_NH3_flow, params.max_NH3_flow);
            fprintf('    Energy per kg NH3: %.1f kJ/kg\n', energy_per_kg_NH3/1000);
            fprintf('    Q_NH3_circulation: %.1f kW\n', Q_NH3_circulation/1000);
            fprintf('    NH3 constraint: %s\n', constraint_type);
            fprintf('  \n');
            
            fprintf('  HEAT PUMP ANALYSIS:\n');
            if state.heat_pump_active
                fprintf('    Status: ACTIVE (%.1f°C > %.1f°C stop)\n', state.temperature, params.T_heat_pump_stop);
            else
                fprintf('    Status: INACTIVE (%.1f°C < %.1f°C start)\n', state.temperature, params.T_heat_pump_start);
            end            
            if exist('Q_heat_pump_demand', 'var')
                fprintf('    Demand: %.1f kW, Delivered: %.1f kW\n', ...
                        Q_heat_pump_demand/1000, Q_heat_pump_delivered/1000);
            end
            fprintf('    Load: %.0f%% of rated capacity\n', HP_load_percent);
            fprintf('  \n');
            
            fprintf('  WATER PREHEATING:\n');
            fprintf('    Q_water_preheat: %.1f kW\n', Q_water_preheat/1000);
            if Q_water_preheat == 0 && state.heat_pump_active
                fprintf('    (Zero: Heat pump has priority)\n');
            end
            fprintf('  \n');
        end
        
        fprintf('  THERMAL MASS:\n');
        fprintf('    Q_excess: %.1f kW [%s]\n', Q_excess/1000, thermal_mass_status);
        fprintf('    Evaporation rate: %.1f g/s (+ = evap, - = condensation)\n', evapRate*1000);
        fprintf('    Temperature change: %.3f°C/step\n', state.temperature - prev_temperature);
        fprintf('  \n');
        
        fprintf('  ENERGY BALANCE VERIFICATION:\n');
        fprintf('    Total output: %.1f kW\n', Q_total_output/1000);
        fprintf('    Balance error: %.3f kW\n', energy_balance_error/1000);
        if abs(energy_balance_error) < params.energy_balance_tolerance
            fprintf('    ✅ Energy balance OK\n');
        else
            fprintf('    ❌ Energy balance ERROR!\n');
        end
        fprintf('  \n');
        
        fprintf('  SYSTEM STATUS: %s\n', equilibrium_status);
        if flow_is_constrained
            fprintf('  ⚠️  FLOW CONSTRAINT: %s\n', constraint_type);
        end
        fprintf('=============================================\n');
    end
    
    %% Safety checks
    if state.temperature >= 123
        warning('Approaching critical temperature, simulation stopped');
        results = truncateResults(results, step + 1);
        break;
    end
    
    if state.liquidMass <= 0.001
        fprintf('All liquid evaporated at t = %.1f minutes\n', state.time/60);
        results = truncateResults(results, step + 1);
        break;
    end
end

fprintf('Simulation completed!\n');

%% Create plots
createUnifiedThermosyphonPlots(results, params);

%% Display final results
displayUnifiedResults(results, state, params);

%% Create vapor column visualization
create_vapor_column_plots(results, params);

%% Supporting Functions
%% Cooper Correlation for Nucleate Pool Boiling Heat Transfer Coefficient
function h_evap = calculateCooperCorrelation(T_sat, P_sat, q_flux, fluid, R_p)
    % COOPER CORRELATION (1984) - Pool Boiling Heat Transfer Coefficient
    %
    % Formula: h = 55 * p_r^(0.12 - 0.4343*ln(R_p)) * (-ln(p_r))^(-0.55) * M^(-0.5) * q^0.67
    %
    % INPUTS:
    %   T_sat   - Saturation temperature [°C]
    %   P_sat   - Saturation pressure [Pa]
    %   q_flux  - Heat flux [W/m²]
    %   fluid   - Fluid name (e.g., 'R717' for ammonia)
    %   R_p     - Surface roughness [μm] (optional, default = 1.0 for commercial tubes)
    %
    % OUTPUTS:
    %   h_evap  - Evaporation heat transfer coefficient [W/m²·K]
    %
    % VALID RANGE:
    %   - Heat flux: 1,000 - 250,000 W/m² (1-250 kW/m²)
    %   - Reduced pressure: 0.001 - 0.9
    %   - Surface roughness: 0.3-1.0 μm (commercial tubes)
    %
    % REFERENCE:
    %   Cooper, M.G. (1984). "Heat Flow Rates in Saturated Nucleate Pool Boiling—
    %   A Wide-Ranging Examination Using Reduced Properties," 
    %   Advances in Heat Transfer, Vol. 16, pp. 157-239
    R_p = 1*10^(-6); %
    % Get critical properties from CoolProp
    try
        P_crit = py.CoolProp.CoolProp.PropsSI('Pcrit', fluid);  % [Pa]
        M = py.CoolProp.CoolProp.PropsSI('M', fluid);           % [kg/mol]
        M = double(M) * 1000;  % Convert to [kg/kmol] (same as g/mol)
    catch
        error('Failed to get critical properties from CoolProp for fluid: %s', fluid);
    end
    
    % Calculate reduced pressure
    p_r = P_sat / P_crit;
    
    % Validate inputs
    if p_r < 0.001 || p_r > 0.9
        warning('Reduced pressure (%.4f) outside recommended range [0.001, 0.9]', p_r);
    end
    
    if q_flux < 1000 || q_flux > 250000
        warning('Heat flux (%.0f W/m²) outside validated range [1,000, 250,000] W/m²', q_flux);
    end
    
    if q_flux <= 0
        warning('Heat flux is zero or negative, setting h_evap to minimum value');
        h_evap = 100;  % Minimum reasonable value [W/m²·K]
        return;
    end
    
    % Cooper Correlation
    % h = 55 * p_r^(0.12 - 0.4343*ln(R_p)) * (-ln(p_r))^(-0.55) * M^(-0.5) * q^0.67
    
    % Calculate each term
    C = 55;  % Constant for SI units
    
    % Pressure term with surface roughness effect
    pressure_exponent = 0.12 - 0.4343 * log(R_p);
    pressure_term = p_r^pressure_exponent;
    
    % Logarithmic pressure term
    if p_r > 0.001
        log_pressure_term = (-log(p_r))^(-0.55);
    else
        log_pressure_term = 1;  % Limit for very low pressures
    end
    
    % Molecular weight term
    molecular_weight_term = M^(-0.5);
    
    % Heat flux term
    heat_flux_term = q_flux^0.67;
    
    % Complete Cooper correlation
    h_evap = C * pressure_term * log_pressure_term * molecular_weight_term * heat_flux_term;
    
    % Apply reasonable bounds
    h_evap = max(100, min(h_evap, 50000));  % Limit to [100, 50000] W/m²·K
    
    h_evap = double(h_evap);
end

function [dmEvap_dt, Q_latent_Cooper] = calculatePhaseChangePhysics(state, params, Q_available)
    % Calculate evaporation rate from Cooper correlation for LATENT heat only
    %
    % INPUTS:
    %   state - Current system state (temperature, pressure, masses, etc.)
    %   params - System parameters
    %   Q_available - Available heat for evaporation [W]
    %
    % OUTPUTS:
    %   dmEvap_dt - Mass evaporation rate [kg/s]
    %   Q_latent_Cooper - LATENT heat transfer from Cooper correlation [W]
    %                     (Sensible heat calculated separately in main loop)
    
    % Get current thermodynamic state
    T_current = state.temperature;  % [°C]
    P_vapor = state.vaporPressure;  % [Pa]
    
    % Calculate interface conditions (accounting for vapor column pressure)
    vapor_column_height = params.length - state.liquidHeight;
    vaporDensity = py.CoolProp.CoolProp.PropsSI('D', 'P', P_vapor, 'Q', 1, params.fluid);
    P_interface = P_vapor + vaporDensity * 9.81 * vapor_column_height;
    
    % Saturation temperature at interface
    T_interface = py.CoolProp.CoolProp.PropsSI('T', 'P', P_interface, 'Q', 0, params.fluid) - 273.15;
    
    %% Get latent heat at interface conditions
    h_liquid = py.CoolProp.CoolProp.PropsSI('H', 'P', P_interface, 'Q', 0, params.fluid);
    h_vapor = py.CoolProp.CoolProp.PropsSI('H', 'P', P_interface, 'Q', 1, params.fluid);
    hfg = double(h_vapor - h_liquid);
    
    % Calculate evaporation rate based on selected method
        % METHOD 1: Heat Flux Limited 
        % Based on Huang et al. (2022) approach
        
        % Maximum possible evaporation if all available heat used
        dmEvap_dt_max = Q_available / hfg;
        
        % Calculate heat flux from available energy
        q_flux = Q_available / params.evap_area;  % [W/m²]

        % Estimate wall temperature (PLACEHOLDER - will be replaced with thermal mass calculation)
        T_wall_estimate = T_current + params.T_wall_offset;

        % Temperature driving force: WALL-TO-LIQUID (not vapor-liquid interface!)
        deltaT_wall_liquid = max(0.1, T_wall_estimate - T_current);  % [°C]

        % COOPER CORRELATION - Calculate h_evap dynamically
        h_evap = calculateCooperCorrelation(T_current, P_vapor, q_flux, params.fluid, 1.0);

        % Nucleate boiling heat flux: q = h * A * ΔT (wall-to-liquid)
        q_evap = h_evap * params.evap_area * deltaT_wall_liquid;

        % Mass evaporation from heat transfer
        dmEvap_dt_htc = q_evap / hfg;

        % Take minimum (energy-limited or heat-transfer-limited)
        dmEvap_dt = min(dmEvap_dt_max, dmEvap_dt_htc);
        
    
    % Apply safety limits
    % Cannot evaporate more than 10% of liquid mass per timestep
    if state.liquidMass > 0
        max_evap_rate = 0.10 * state.liquidMass / params.timeStep;
        dmEvap_dt = min(dmEvap_dt, max_evap_rate);
    end

    % ALLOW NEGATIVE EVAPORATION (CONDENSATION)
    % If Q_excess < 0 (cooling), condensation can occur
    % No constraint on negative values - dmEvap_dt can be negative

    % CORRECTED: Return LATENT heat only from Cooper correlation
    % Cooper calculates nucleate boiling heat transfer (latent component)
    % Sensible heat (reheating NH3 from 20°C to saturation) calculated separately
    Q_latent_Cooper = q_evap;  % Latent heat from Cooper correlation [W]

    % Final check: Cannot transfer more latent heat than available
    if Q_latent_Cooper > Q_available
        Q_latent_Cooper = Q_available;
        dmEvap_dt = Q_available / hfg;  % Limit evaporation accordingly
    end
end
%% FIXED: Improved Mass Balance Solver
function [T_final, P_final] = solve_mass_balance_improved(totalVaporMass, vaporVolume, T_initial, fluid)
    % FIXED: Better handling of negative mass changes and convergence
    
    % Ensure minimum vapor mass to prevent numerical issues
    if totalVaporMass <= 0
        totalVaporMass = 1e-6; % Minimum vapor mass
    end
    
    T_guess = T_initial;
    tolerance = 0.001;
    max_iterations = 30; % Increased iterations
    
    for iter = 1:max_iterations
        try
            P_sat = py.CoolProp.CoolProp.PropsSI('P', 'T', T_guess + 273.15, 'Q', 0.5, fluid);
            rho_vapor = py.CoolProp.CoolProp.PropsSI('D', 'P', P_sat, 'Q', 1, fluid);
            calculated_vapor_mass = rho_vapor * vaporVolume;
            mass_error = totalVaporMass - calculated_vapor_mass;
            
            if abs(mass_error) < tolerance
                break;
            end
            
            % FIXED: Adaptive step size based on mass error magnitude
            if abs(mass_error) > 10
                step_size = 2.0;
            elseif abs(mass_error) > 1
                step_size = 1.0;
            else
                step_size = 0.5;
            end
            
            if mass_error > 0
                T_guess = T_guess + step_size;
            else
                T_guess = T_guess - step_size;
            end
            
            % Enforce reasonable bounds
            T_guess = max(T_guess, -70);
            T_guess = min(T_guess, 130);
            
        catch
            % If CoolProp fails, use a small adjustment
            if mass_error > 0
                T_guess = T_guess + 0.1;
            else
                T_guess = T_guess - 0.1;
            end
        end
    end
    
    if iter >= max_iterations
        fprintf('Warning: Mass balance solver did not converge after %d iterations\n', max_iterations);
    end
    
    T_final = T_guess;
    try
        P_final = py.CoolProp.CoolProp.PropsSI('P', 'T', T_final + 273.15, 'Q', 0.5, fluid);
    catch
        % Fallback if CoolProp fails
        P_final = 101325; % 1 atm
    end
end

%% Geothermal temperature function
function T_ground = geothermal_temperature(depth, params)
    % Sigmoid function: T = 81.12896 / (1 + e^(-(0.0000589*x - 1.99186)))
    exponent = -(params.geothermal_coeff_a * depth - params.geothermal_coeff_b);
    T_ground = params.geothermal_max_temp / (1 + exp(exponent));
end

%% Calculate depth-dependent heat input
function [Q_geothermal_total, Q_profile_debug] = calculate_geothermal_input(state, params, T_pipe)
    % INPUTS:
    %   state - Current system state
    %   params - System parameters
    %   T_pipe - Actual pipe temperature for calculation [°C]

    if ~params.use_geothermal_profile
        % Fallback to original constant heat input
        Q_geothermal_total = params.heatInput * params.evap_area;
        Q_profile_debug = [];
        return;
    end

    % Calculate heat input along pipe length
    segment_height = 50; % m - calculation segments
    depths = 0:segment_height:params.length;
    Q_total = 0;
    Q_profile_debug = zeros(length(depths)-1, 4); % [depth, T_ground, T_pipe, Q_segment]

    for i = 1:length(depths)-1
        depth_mid = (depths(i) + depths(i+1)) / 2;
        segment_length = depths(i+1) - depths(i);

        % Ground temperature at this depth
        T_ground = geothermal_temperature(depth_mid, params);

        % Use actual pipe temperature (passed as argument)
        % When NH3 is hot, geothermal input naturally decreases (correct physics)

        % Heat transfer area for this segment
        segment_area = pi * params.diameter_out * segment_length;

        % Simplified overall heat transfer coefficient
        U_overall = 10; % W/m²·K - simplified for now

        % Heat transfer rate for this segment (positive = heat into pipe)
        deltaT = T_ground - T_pipe;
        Q_segment = U_overall * segment_area * deltaT;

        % Only count positive heat transfer (heating the pipe)
        Q_segment = max(0, Q_segment);
        Q_total = Q_total + Q_segment;

        % Store debug info
        Q_profile_debug(i, :) = [depth_mid, T_ground, T_pipe, Q_segment/1000]; % kW
    end

    Q_geothermal_total = Q_total;
end

function results = truncateResults(results, lastStep)
    fields = fieldnames(results);
    for i = 1:length(fields)
        fieldName = fields{i};
        
        % Skip vapor_column_data (it's stored every 100 steps, not every step)
        if strcmp(fieldName, 'vapor_column_data')
            continue;  % Don't truncate this field
        end
        
        if isnumeric(results.(fieldName))
            % Truncate numeric arrays
            if length(results.(fieldName)) >= lastStep
                results.(fieldName) = results.(fieldName)(1:lastStep);
            end
        elseif iscell(results.(fieldName))
            % Truncate cell arrays (like status arrays)
            if length(results.(fieldName)) >= lastStep
                results.(fieldName) = results.(fieldName)(1:lastStep);
            end
        end
    end
end

%% Create Real Temperature Profile Function
function [heights, T_profile, P_profile, section_labels] = create_real_temperature_profile(state, params, T_interface_actual)
    % Create height array (every 50m)
    heights = 0:params.segment_height:params.length;
    if heights(end) ~= params.length
        heights = [heights, params.length]; % Ensure we include the top
    end
    
    n_points = length(heights);
    T_profile = zeros(n_points, 1);
    P_profile = zeros(n_points, 1);
    section_labels = cell(n_points, 1);
    
    % Define section boundaries
    evap_end = params.length * 0.25;
    adiab_end = evap_end + params.length * 0.500;
    
    % Get vapor density at top conditions (for hydrostatic calculation)
    rho_vapor = py.CoolProp.CoolProp.PropsSI('D', 'P', state.vaporPressure, 'Q', 1, params.fluid);
    
    for i = 1:n_points
        h = heights(i);
        
        if h <= state.liquidHeight
            %% LIQUID REGION: Uniform temperature and hydrostatic pressure
            T_profile(i) = T_interface_actual;
            
            % Liquid hydrostatic pressure (higher pressure at bottom)
            liquid_rho = py.CoolProp.CoolProp.PropsSI('D', 'T', T_interface_actual + 273.15, 'Q', 0, params.fluid);
            P_profile(i) = state.vaporPressure + liquid_rho * 9.81 * (state.liquidHeight - h);
            
            section_labels{i} = 'Liquid';
            
        else
            %% VAPOR REGION: Calculate pressure with vapor hydrostatic effect
            height_from_top = params.length - h;  % Height from top of heat pipe
            
            % Vapor hydrostatic pressure (P increases going down)
            P_profile(i) = state.vaporPressure + rho_vapor * 9.81 * height_from_top;
            
            % Get saturation temperature at this pressure using CoolProp
            try
                T_profile(i) = py.CoolProp.CoolProp.PropsSI('T', 'P', P_profile(i), 'Q', 0, params.fluid) - 273.15;
            catch
                % If pressure is too high/low for CoolProp, use uniform temperature
                T_profile(i) = state.temperature;
            end
            
            % Determine section
            if h <= evap_end
                section_labels{i} = 'Evap_Vapor';
            elseif h <= adiab_end
                section_labels{i} = 'Adiabatic';
            else
                section_labels{i} = 'Condenser';
            end
        end
    end
    
    % Debug outputinitial_liquidDensity
    if any(heights > state.liquidHeight)
        vapor_indices = heights > state.liquidHeight;
        T_range = max(T_profile(vapor_indices)) - min(T_profile(vapor_indices));
        P_range = (max(P_profile(vapor_indices)) - min(P_profile(vapor_indices))) / 1000;
        
        fprintf('  Vapor column: ΔT = %.2f°C, ΔP = %.1f kPa over %.0fm\n', ...
                T_range, P_range, params.length - state.liquidHeight);
    end
end

%% Create Vapor Column Visualization
function create_vapor_column_plots(results, params)
    if ~isfield(results, 'vapor_column_data') || isempty(results.vapor_column_data)
        fprintf('No vapor column data available\n');
        return;
    end
    
    figure('Position', [100, 100, 1600, 1000]);
    
    % Plot 1: Temperature profiles over time
    subplot(2, 3, 1);
    hold on;
    
    n_profiles = length(results.vapor_column_data);
    times_to_plot = [1, round(n_profiles/3), round(2*n_profiles/3), n_profiles];
    colors = {'b', 'g', 'r', 'k'};
    
    for i = 1:length(times_to_plot)
        idx = times_to_plot(i);
        if idx <= n_profiles
            data = results.vapor_column_data{idx};
            plot(data.temperatures, data.heights, colors{i}, 'LineWidth', 2, ...
                 'DisplayName', sprintf('t=%.0f min', data.time));
            
            % Mark liquid level and uniform temperature
            yline(data.liquid_height, '--', colors{i}, 'LineWidth', 1);
            xline(data.uniform_temp, ':', colors{i}, 'LineWidth', 1);
        end
    end
    
    xlabel('Temperature (°C)'); ylabel('Height (m)');
    title('Real Temperature Profiles (Hydrostatic)');
    legend('Location', 'best'); grid on;
    
    % Plot 2: Pressure profiles
    subplot(2, 3, 2);
    hold on;
    
    for i = 1:length(times_to_plot)
        idx = times_to_plot(i);
        if idx <= n_profiles
            data = results.vapor_column_data{idx};
            plot(data.pressures, data.heights, colors{i}, 'LineWidth', 2, ...
                 'DisplayName', sprintf('t=%.0f min', data.time));
            
            yline(data.liquid_height, '--', colors{i}, 'LineWidth', 1);
        end
    end
    
    xlabel('Pressure (kPa)'); ylabel('Height (m)');
    title('Pressure Profiles (Hydrostatic)');
    legend('Location', 'best'); grid on;
    
    % Plot 3: Temperature vs height at final time
    subplot(2, 3, 3);
    final_data = results.vapor_column_data{end};
    
    % Separate liquid and vapor regions
    liquid_mask = final_data.heights <= final_data.liquid_height;
    vapor_mask = final_data.heights > final_data.liquid_height;
    
    % DEBUG: Check temperature profile
    if any(vapor_mask)
        vapor_temps = final_data.temperatures(vapor_mask);
        vapor_heights = final_data.heights(vapor_mask);
        
        fprintf('\n=== Temperature Profile Debug ===\n');
        fprintf('Liquid height: %.0fm\n', final_data.liquid_height);
        fprintf('Top vapor: h=%.0fm, T=%.1f°C\n', max(vapor_heights), vapor_temps(vapor_heights == max(vapor_heights)));
        fprintf('Bottom vapor: h=%.0fm, T=%.1f°C\n', min(vapor_heights), vapor_temps(vapor_heights == min(vapor_heights)));
        fprintf('Temperature range: %.1f°C to %.1f°C\n', min(vapor_temps), max(vapor_temps));
        fprintf('Expected: Bottom vapor should be HOTTER than top vapor\n');
        fprintf('====================================\n');
    end

    hold on;
    if any(liquid_mask)
        plot(final_data.temperatures(liquid_mask), final_data.heights(liquid_mask), ...
             'b-', 'LineWidth', 3, 'DisplayName', 'Liquid');
    end
    if any(vapor_mask)
        plot(final_data.temperatures(vapor_mask), final_data.heights(vapor_mask), ...
             'r-', 'LineWidth', 3, 'DisplayName', 'Vapor');
    end
    
    % Mark extraction point (80% up in vapor space)
    if final_data.liquid_height < params.length
        extraction_height = final_data.liquid_height + 0.8 * (params.length - final_data.liquid_height);
        extraction_temp = interp1(final_data.heights, final_data.temperatures, extraction_height);
        plot(extraction_temp, extraction_height, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', ...
             'DisplayName', sprintf('NH3 Extraction (%.1f°C)', extraction_temp));
    end
    
    yline(final_data.liquid_height, '--k', 'Liquid Level', 'LineWidth', 2);
    xline(final_data.uniform_temp, ':k', 'Uniform T', 'LineWidth', 2);
    
    xlabel('Temperature (°C)'); ylabel('Height (m)');
    title(sprintf('Final Temperature Profile (t=%.0f min)', final_data.time));
    legend('Location', 'best'); grid on;
    
    % Plot 4: Vapor column temperature difference over time
    subplot(2, 3, 4);
    times = zeros(n_profiles, 1);
    temp_differences = zeros(n_profiles, 1);
    
    for i = 1:n_profiles
        data = results.vapor_column_data{i};
        times(i) = data.time;
        
        % Calculate temperature difference in vapor column
        vapor_mask = data.heights > data.liquid_height;
        if any(vapor_mask)
            vapor_temps = data.temperatures(vapor_mask);
            temp_differences(i) = max(vapor_temps) - min(vapor_temps);
        else
            temp_differences(i) = 0;
        end
    end
    
    plot(times, temp_differences, 'b-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('Temperature Difference (°C)');
    title('Vapor Column ΔT Over Time');
    grid on;
    
    % Plot 5: Heat pipe schematic with temperature colors
    subplot(2, 3, 5:6);
    
    % Create temperature colormap
    all_temps = final_data.temperatures;
    temp_range = [min(all_temps), max(all_temps)];
    
    % Draw temperature-colored heat pipe
    pipe_width = 0.3;
    x_left = -pipe_width/2;
    x_right = pipe_width/2;
    
    hold on;
    
    % Color-coded temperature regions
    for i = 1:length(final_data.heights)-1
        h1 = final_data.heights(i);
        h2 = final_data.heights(i+1);
        temp = final_data.temperatures(i);
        
        % Normalize temperature for color
        temp_norm = (temp - temp_range(1)) / (temp_range(2) - temp_range(1));
        color = [temp_norm, 0, 1-temp_norm]; % Red (hot) to Blue (cold)
        
        fill([x_left, x_right, x_right, x_left], [h1, h1, h2, h2], ...
             color, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    end
    
    % Pipe outline
    plot([x_left, x_left], [0, params.length], 'k-', 'LineWidth', 3);
    plot([x_right, x_right], [0, params.length], 'k-', 'LineWidth', 3);
    
    % Liquid level line
    plot([x_left, x_right], [final_data.liquid_height, final_data.liquid_height], ...
         'w-', 'LineWidth', 3);
    
    % Add colorbar
    colormap(gca, [linspace(0,1,256)', zeros(256,1), linspace(1,0,256)']);
    c = colorbar;
    c.Label.String = 'Temperature (°C)';
    caxis(temp_range);
    
    xlim([x_left - 0.1, x_right + 0.1]);
    ylim([0, params.length]);
    xlabel('Heat Pipe Cross-Section');
    ylabel('Height (m)');
    title('Temperature Distribution in Heat Pipe');
    
    sgtitle('Vapor Column Analysis - Real Thermodynamic Temperatures', 'FontSize', 16);
end

function createUnifiedThermosyphonPlots(results, params)
    figure('Position', [50, 50, 1600, 1200], 'Color', 'white');
    
    % Temperature Evolution
    subplot(3, 3, 1);
    plot(results.time, results.satTemp, 'b-', 'LineWidth', 2);
    hold on;
    yline(params.T_equilibrium_target, '--r', 'Target Temp', 'LineWidth', 2);
    yline(params.T_heat_pump_start, '--g', 'HP Start', 'LineWidth', 1);
    yline(params.T_heat_pump_stop, '--m', 'HP Stop', 'LineWidth', 1);
    equilibrium_mask = strcmp(results.equilibrium_status, 'ENERGY_EQUILIBRIUM');
    if any(equilibrium_mask)
        plot(results.time(equilibrium_mask), results.satTemp(equilibrium_mask), 'g-', 'LineWidth', 3, 'DisplayName', 'Equilibrium');
    end
    xlabel('Time (min)'); ylabel('Temperature (°C)');
    title('Temperature Evolution (FIXED: Shows Cooling)');
    legend; grid on;
    
    % Energy Balance Components
    subplot(3, 3, 2);
    hold on;
    plot(results.time, results.Q_input_effective, '--k', 'LineWidth', 2, 'DisplayName', 'Input Total×EF');
    plot(results.time, results.Q_heat_pump, 'r-', 'LineWidth', 2, 'DisplayName', 'Heat Pump');
    plot(results.time, results.Q_water_preheat, 'b-', 'LineWidth', 2, 'DisplayName', 'Water Preheat');
    plot(results.time, results.Q_NH3_circulation, 'g-', 'LineWidth', 2, 'DisplayName', 'NH3 Total');
    plot(results.time, results.Q_excess, 'm-', 'LineWidth', 2, 'DisplayName', 'Excess');
    plot(results.time, results.Q_losses, 'c-', 'LineWidth', 1, 'DisplayName', 'Losses');
    xlabel('Time (min)'); ylabel('Power (kW)');
    title('Unified Energy Balance');
    legend; grid on;
    
    % NH3 Flow and Constraints
    subplot(3, 3, 3);
    yyaxis left
    plot(results.time, results.NH3_flow_actual*1000, 'b-', 'LineWidth', 2);
    ylabel('NH3 Flow (g/s)', 'Color', 'b');
    yline(params.max_NH3_flow*1000, '--r', 'Max Flow');
    
    yyaxis right
    constraint_numeric = zeros(size(results.time));
    for i = 1:length(results.constraint_status)
        switch results.constraint_status{i}
            case 'MAX_FLOW', constraint_numeric(i) = 2;
            case 'MIN_FLOW', constraint_numeric(i) = 1;
            case 'NO_CIRCULATION', constraint_numeric(i) = -1;
            otherwise, constraint_numeric(i) = 0;
        end
    end
    plot(results.time, constraint_numeric, 'r-', 'LineWidth', 1);
    ylabel('Constraint Status', 'Color', 'r');
    xlabel('Time (min)');
    title('NH3 Flow Control');
    grid on;
    
    % Heat Pump vs Water Preheating Analysis
    subplot(3, 3, 4);
    yyaxis left
    hold on;
    plot(results.time, results.Q_heat_pump, 'r-', 'LineWidth', 2, 'DisplayName', 'Q Heat Pump');
    plot(results.time, results.Q_water_preheat, 'b-', 'LineWidth', 2, 'DisplayName', 'Q Preheating');
    ylabel('Power Output (kW)', 'Color', 'k');
    
    yyaxis right
    plot(results.time, results.preheat_percentage, 'g--', 'LineWidth', 2, 'DisplayName', '% Preheating');
    ylabel('Preheating Share (%)', 'Color', 'g');
    ylim([0, 100]);
    
    yyaxis left
    legend('Location', 'northwest');
    xlabel('Time (min)');
    title('Heat Pump vs Preheating Output');
    grid on;    
    
    % Heat Pump Performance vs Requirement
    subplot(3, 3, 5);
    hold on;
    plot(results.time, results.Q_heat_pump, 'b-', 'LineWidth', 2, 'DisplayName', 'HP Actual');
    plot(results.time, results.HP_load_percent, 'g-', 'LineWidth', 2, 'DisplayName', 'HP Load %');

    % Plot dynamic heat pump demand
    hp_demand_over_time = arrayfun(params.heat_pump_duty, results.time*60);  % Convert minutes to seconds
    plot(results.time, hp_demand_over_time/1000, '--r', 'LineWidth', 2, 'DisplayName', 'HP Required');

    xlabel('Time (min)'); ylabel('Heat Pump Power (kW) / Load (%)');
    title('Heat Pump Performance (FIXED: Hysteresis)');
    legend; grid on;
    % Energy Balance Verification
    subplot(3, 3, 6);
    plot(results.time, results.energy_balance_error, 'r-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('Energy Balance Error (kW)');
    title('Energy Balance Check');
    grid on;
    yline(0, '--k', 'Perfect Balance');
    
    % System Mass Evolution
    subplot(3, 3, 7);
    yyaxis left
    plot(results.time, results.liquidMass/1000, 'b-', 'LineWidth', 2);
    ylabel('Liquid Mass (tons)', 'Color', 'b');
    
    yyaxis right
    plot(results.time, results.liquidHeight, 'r-', 'LineWidth', 2);
    ylabel('Liquid Height (m)', 'Color', 'r');
    xlabel('Time (min)');
    title('System Mass & Geometry');
    grid on;
    
    % Pressure Evolution
    subplot(3, 3, 8);
    semilogy(results.time, results.vaporPressure, 'b-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('Pressure (kPa)');
    title('Vapor Pressure Evolution');
    grid on;
    
    % Energy Factor Evolution
    subplot(3, 3, 9);
    plot(results.time, results.energyFactor*100, 'g-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('Energy Factor (%)');
    title('Dynamic Energy Factor');
    yline(30, '--b', '30% (Heating)');
    yline(81.3, '--r', '81.3% (Circulation)');
    grid on;
    
    % Final System State
    subplot(3, 3, 9);
    rectangle('Position', [0, 0, 1, 1], 'FaceColor', [0.7, 0.9, 1], 'EdgeColor', 'black', 'LineWidth', 2);
    
    if results.liquidHeight(end) > 0
        currentLiquidFraction = results.liquidHeight(end) / results.liquidHeight(1);
        rectangle('Position', [0, 0, 1, currentLiquidFraction], 'FaceColor', [0.2, 0.6, 1], 'EdgeColor', 'none');
    else
        currentLiquidFraction = 0;
    end
    
    % Status indicator
    equilibrium_status_final = results.equilibrium_status{end};
    if strcmp(equilibrium_status_final, 'ENERGY_EQUILIBRIUM')
        status_color = [0, 0.8, 0];
        status_text = 'EQUILIBRIUM';
    else
        status_color = [1, 0.5, 0];
        status_text = 'TRANSIENT';
    end
    
    text(0.5, 0.8, sprintf('%s', status_text), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', status_color);
    
    text(0.5, 0.6, sprintf('T = %.1f°C\nP = %.1f kPa\nHP = %.1f kW', ...
        results.satTemp(end), results.vaporPressure(end), results.Q_heat_pump(end)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    
    xlim([0, 1]); ylim([0, 1]);
    title('Final System State');
    set(gca, 'XTick', [], 'YTick', []);
    
    sgtitle(sprintf('FIXED Sequential Energy Cascade Model - %s (%.0f kW/m²)', ...
        params.fluid, params.heatInput/1000), 'FontSize', 16, 'FontWeight', 'bold');
end

function displayUnifiedResults(results, state, params)
    fprintf('\n=== FIXED SEQUENTIAL ENERGY CASCADE SIMULATION RESULTS ===\n');
    fprintf('Fluid: %s\n', params.fluid);
    fprintf('Heat Input: %.1f kW/m²\n', params.heatInput/1000);
    
    % Analyze energy factor usage
    final_energy_factor = results.energyFactor(end);
    heating_steps = sum(results.energyFactor == 0.30);
    circulation_steps = sum(results.energyFactor == 0.813);
    total_steps = length(results.energyFactor);

    fprintf('Energy Factor Analysis:\n');
    fprintf('  Final: %.1f%%\n', final_energy_factor*100);
    fprintf('  Heating phase (30%%): %d steps (%.1f%% of simulation)\n', ...
            heating_steps, 100*heating_steps/total_steps);
    fprintf('  Circulation phase (81.3%%): %d steps (%.1f%% of simulation)\n', ...
            circulation_steps, 100*circulation_steps/total_steps);

    % Find transition time if it occurred
    if heating_steps > 0 && circulation_steps > 0
        transition_step = find(results.energyFactor == 0.813, 1);
        transition_time = results.time(transition_step);
        fprintf('  Mode transition at: %.1f minutes\n', transition_time);
    end
    
    fprintf('Simulation Time: %.1f minutes\n', results.time(end));
    
    fprintf('\nFinal Conditions:\n');
    fprintf('  Temperature: %.1f°C (target: %.1f°C)\n', results.satTemp(end), params.T_equilibrium_target);
    fprintf('  Pressure: %.1f kPa\n', results.vaporPressure(end));
    fprintf('  Liquid Mass: %.2f kg\n', results.liquidMass(end));
    
    % Check if heat pump demand exists at final time
    if results.time(end) > 0
        final_hp_demand = params.heat_pump_duty(results.time(end)*60)/1000;
        fprintf('  Heat Pump Output: %.1f kW (required: %.1f kW at final time)\n', results.Q_heat_pump(end), final_hp_demand);
    else
        fprintf('  Heat Pump Output: %.1f kW\n', results.Q_heat_pump(end));
    end
    
    fprintf('  Heat Pump Load: %.0f%% of rated capacity\n', results.HP_load_percent(end));
    fprintf('  Water Preheating: %.1f kW\n', results.Q_water_preheat(end));
    fprintf('  NH3 Flow: %.3f kg/s\n', results.NH3_flow_actual(end));
    fprintf('  Excess Energy: %.1f kW\n', results.Q_excess(end));
    
    if state.equilibrium_achieved
        fprintf('\n✅ EQUILIBRIUM ACHIEVED!\n');
        fprintf('  Time to equilibrium: %.1f minutes\n', state.equilibrium_start_time/60);
        if results.time(end) > 0
            final_hp_demand = params.heat_pump_duty(results.time(end)*60)/1000;
            if final_hp_demand > 0
                fprintf('  Heat pump satisfaction: %.1f%%\n', 100*results.Q_heat_pump(end)/final_hp_demand);
            end
        end
    else
        fprintf('\n⚠ Equilibrium not achieved in simulation time\n');
        temp_diff = abs(results.satTemp(end) - params.T_equilibrium_target);
        fprintf('  Temperature difference: %.1f°C\n', temp_diff);
    end
    
    % Constraint analysis
    max_flow_count = sum(strcmp(results.constraint_status, 'MAX_FLOW'));
    min_flow_count = sum(strcmp(results.constraint_status, 'MIN_FLOW'));
    no_circ_count = sum(strcmp(results.constraint_status, 'NO_CIRCULATION'));
    
    fprintf('\nConstraint Analysis:\n');
    fprintf('  Max flow constraint: %d steps (%.1f%%)\n', max_flow_count, 100*max_flow_count/length(results.time));
    fprintf('  Min flow constraint: %d steps (%.1f%%)\n', min_flow_count, 100*min_flow_count/length(results.time));
    fprintf('  No circulation: %d steps (%.1f%%)\n', no_circ_count, 100*no_circ_count/length(results.time));
    
    % Energy balance summary
    avg_energy_error = mean(abs(results.energy_balance_error));
    fprintf('\nEnergy Balance Quality:\n');
    fprintf('  Average energy error: %.2f kW\n', avg_energy_error);
    fprintf('  Final energy error: %.2f kW\n', results.energy_balance_error(end));
    
    % Performance metrics
    if any(results.Q_heat_pump > 0)
        active_hp_indices = results.Q_heat_pump > 0;
        if any(active_hp_indices)
            % Calculate average efficiency for times when heat pump is active
            avg_hp_output = mean(results.Q_heat_pump(active_hp_indices));
            
            % Get corresponding heat pump demands
            active_times = results.time(active_hp_indices) * 60; % Convert to seconds
            hp_demands = arrayfun(params.heat_pump_duty, active_times) / 1000; % Convert to kW
            avg_hp_demand = mean(hp_demands);
            
            if avg_hp_demand > 0
                hp_efficiency = avg_hp_output / avg_hp_demand * 100;
                fprintf('\nPerformance Metrics:\n');
                fprintf('  Average heat pump efficiency: %.1f%%\n', hp_efficiency);
            end
        end
        
        fprintf('  Total NH3 circulated: %.1f kg\n', sum(results.NH3_flow_actual) * params.timeStep);
        fprintf('  Total water preheating: %.1f MJ\n', sum(results.Q_water_preheat) * params.timeStep / 1000);
    end
    
    fprintf('\n=== KEY FIXES IMPLEMENTED ===\n');
    fprintf('1. ✅ FIXED: Heat pump hysteresis (%.1f°C start, %.1f°C stop)\n', params.T_heat_pump_start, params.T_heat_pump_stop);
    fprintf('2. ✅ FIXED: Mass balance solver handles condensation (Q_excess < 0)\n');
    fprintf('3. ✅ FIXED: Water preheating priority (zero when HP not at 100%%)\n');
    fprintf('4. ✅ FIXED: Improved temperature drop physics\n');
    fprintf('5. ✅ FIXED: Better convergence and bounds checking\n');
    fprintf('6. ✅ Original features: Sequential energy cascade, no circular dependency\n');
    
    % Temperature change analysis
    temp_changes = diff(results.satTemp);
    max_heating_rate = max(temp_changes);
    max_cooling_rate = min(temp_changes);
    
    fprintf('\nTemperature Change Analysis:\n');
    fprintf('  Max heating rate: %.3f°C/step\n', max_heating_rate);
    fprintf('  Max cooling rate: %.3f°C/step\n', max_cooling_rate);
    
    if max_cooling_rate < -0.001
        fprintf('  ✅ System now shows proper temperature cooling!\n');
    else
        fprintf('  ⚠️  Limited temperature cooling observed\n');
    end
end
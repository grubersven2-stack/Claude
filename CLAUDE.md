# CLAUDE.md - Thermosyphon Simulation Codebase

## Repository Overview

This repository contains a sophisticated MATLAB-based simulation of a **Unified Thermosyphon Model with Continuous NH3 Circulation**. The model simulates a geothermal heat pipe system coupled with a heat pump and water preheating functionality.

### Purpose

The simulation models the thermodynamic behavior of a deep geothermal heat pipe (1800m) using ammonia (NH3/R717) as the working fluid. It includes:
- Two-phase flow dynamics (liquid and vapor)
- Heat pump coupling with hysteresis control
- Water preheating capabilities
- Energy balance verification
- Real-time visualization of system performance

## Codebase Structure

### File Organization

```
Claude/
├── Energy1.m           # Main simulation script (1381 lines)
└── CLAUDE.md          # This file
```

### Code Architecture

The `Energy1.m` file is organized into the following sections:

1. **Initialization (Lines 1-195)**: System parameters, CoolProp setup, initial conditions
2. **Main Simulation Loop (Lines 196-623)**: Time-stepping solver with energy cascade
3. **Supporting Functions (Lines 636-1381)**: Physics models and visualization

#### Key Functions

| Function | Line Range | Purpose |
|----------|------------|---------|
| `calculatePhaseChangePhysics()` | 638-712 | Calculate evaporation rate from heat transfer physics |
| `solve_mass_balance_improved()` | 714-777 | Iterative solver for thermodynamic equilibrium |
| `geothermal_temperature()` | 780-784 | Sigmoid geothermal temperature profile |
| `calculate_geothermal_input()` | 787-830 | Depth-dependent heat input calculation |
| `truncateResults()` | 832-854 | Truncate results arrays if simulation stops early |
| `create_real_temperature_profile()` | 857-924 | Generate height-dependent temperature profile |
| `create_vapor_column_plots()` | 927-1099 | Visualize vapor column behavior |
| `createUnifiedThermosyphonPlots()` | 1101-1258 | Generate comprehensive system plots |
| `displayUnifiedResults()` | 1260-1381 | Print detailed simulation results |

## Physical Model

### System Description

The model simulates a **1800m vertical heat pipe** divided into three sections:
- **Evaporator** (25% of length): Bottom section where geothermal heat vaporizes NH3
- **Adiabatic section** (50% of length): Middle transport section
- **Condenser** (25% of length): Top section where vapor condenses, releasing heat

### Dimensions
- Total length: 1800 m
- Inner diameter: 0.1778 m
- Outer diameter: 0.1900 m
- Wall material: Steel (7850 kg/m³, 480 J/kg·K)

### Working Fluid: Ammonia (R717)

Ammonia properties are calculated dynamically using **CoolProp** Python library:
- Saturation properties (temperature, pressure)
- Liquid/vapor densities
- Enthalpies and latent heat (hfg)
- Phase change behavior

### Energy Flow Model

The simulation implements a **Sequential Energy Cascade** (no circular dependencies):

```
Q_input_total (Geothermal)
    ↓
Q_input_net = Q_input_total × Energy_Factor - Q_losses
    ↓
├─→ Q_NH3_circulation (NH3 flow × energy_per_kg)
│       ↓
│   ├─→ Q_heat_pump (Priority load)
│   └─→ Q_water_preheat (Excess after HP satisfied)
│
└─→ Q_excess (Thermal mass heating/cooling)
        ↓
    Evaporation or Condensation
```

### Key Physics

1. **Geothermal Temperature Profile** (Line 780)
   ```
   T_ground(depth) = 81.13 / (1 + e^(-(0.0060589*depth - 1.99186)))
   ```

2. **Hydrostatic Pressure Effects**
   - Vapor column pressure: P_interface = P_vapor + ρ_vapor × g × h_vapor
   - Liquid column pressure: P_liquid = P_interface + ρ_liquid × g × h_liquid

3. **Phase Change**
   - Heat flux limited: `q = h × A × ΔT`
   - Evaporation rate: `dm/dt = q / hfg`

4. **Mass Balance Solver** (Line 714)
   - Iterative Newton-like method
   - Ensures total vapor mass equals ρ_vapor × V_vapor
   - Adaptive step size based on error magnitude

## Operating Modes

### 1. Heating Phase (T < T_target)
- Heat pump inactive
- All input energy heats thermal mass
- No NH3 circulation
- Energy factor: 30%

### 2. Circulation Phase (T ≥ T_target OR heat pump active)
- NH3 circulation active (min: 0.01 kg/s, max: 0.2 kg/s)
- Heat pump operates with hysteresis:
  - **Start**: T ≥ 59°C
  - **Stop**: T ≤ 55°C
- Energy factor: 81.3%
- Water preheating gets excess after heat pump satisfied

## Key Parameters

### System Configuration (Lines 16-98)

| Parameter | Variable | Value | Unit | Description |
|-----------|----------|-------|------|-------------|
| Heat input | `heatInput` | 1000 | W/m² | Geothermal heat flux |
| Initial temp | `initialTemp` | 0 | °C | Starting temperature |
| Pipe length | `length` | 1800 | m | Total heat pipe height |
| Liquid height | `initialLiquidHeight` | 360 | m | Initial liquid level (20%) |
| Time step | `timeStep` | 30 | s | Simulation timestep |
| Total time | `totalTime` | 86400 | s | Simulation duration (24h) |

### Heat Pump Parameters (Lines 81-88)

| Parameter | Variable | Value | Unit |
|-----------|----------|-------|------|
| Heat pump demand | `heat_pump_duty` | 40 kW + 30 kW × sin(π×t/24h) | W |
| Target temperature | `T_equilibrium_target` | 60 | °C |
| Return temperature | `T_return_target` | 20 | °C |
| Start threshold | `T_heat_pump_start` | 59 | °C |
| Stop threshold | `T_heat_pump_stop` | 55 | °C |

### NH3 Flow Control (Lines 90-93)

| Parameter | Variable | Value | Unit |
|-----------|----------|-------|------|
| Minimum flow | `min_NH3_flow` | 0.01 | kg/s |
| Maximum flow | `max_NH3_flow` | 0.2 | kg/s |
| Temperature tolerance | `equilibrium_tolerance` | 1.0 | °C |

## Development Workflow

### Prerequisites

1. **MATLAB** (tested version unknown, likely R2020a or later)
2. **Python with CoolProp**:
   ```bash
   pip install CoolProp
   ```
3. MATLAB configured to use Python:
   ```matlab
   pyversion  % Check Python configuration
   ```

### Running the Simulation

```matlab
% In MATLAB command window or terminal:
cd /path/to/Claude
Energy1
```

The script will:
1. Initialize CoolProp
2. Display wall thermal mass parameters
3. Run 24-hour simulation with progress indicators
4. Generate 9-subplot figure with comprehensive results
5. Create vapor column visualization
6. Print detailed results summary

### Typical Execution Time

- **Time steps**: 2880 (24 hours at 30s intervals)
- **Execution time**: ~2-5 minutes (depends on system)
- **Progress updates**: Every 5% (20 milestones)

### Output

The simulation produces:
1. **Figure 1**: 9-panel system analysis
   - Temperature evolution
   - Energy balance components
   - NH3 flow control
   - Heat pump vs preheating
   - Heat pump performance
   - Energy balance error
   - Mass and geometry
   - Pressure evolution
   - Final system state

2. **Figure 2**: Vapor column visualization
   - Temperature profiles over time
   - Pressure profiles
   - Final temperature distribution
   - Temperature differences
   - Heat pipe schematic with colors

3. **Console output**: Detailed results and diagnostics

## Code Conventions

### Commenting Style

- `%%` - Major section headers
- `%` - Inline comments
- `%{` ... `%}` - Multi-line comment blocks (some sections commented out)

### Naming Conventions

| Convention | Example | Usage |
|------------|---------|-------|
| `camelCase` | `heatInput`, `vaporPressure` | Variables and parameters |
| `snake_case` | `heat_pump_duty`, `T_equilibrium_target` | Physical parameters |
| `UPPER_CASE` | `HEATING`, `COOLING`, `ENERGY_EQUILIBRIUM` | Status strings |
| Prefix `Q_` | `Q_input_net`, `Q_heat_pump` | Energy/power quantities |
| Prefix `T_` | `T_interface`, `T_ground` | Temperature variables |
| Prefix `P_` | `P_interface`, `P_vapor` | Pressure variables |
| Prefix `m_` | `m_water_preheat` | Mass or mass flow rates |

### Units Convention

**All internal calculations use SI units:**
- Length: meters [m]
- Mass: kilograms [kg]
- Time: seconds [s]
- Temperature: Celsius [°C] (converted to Kelvin for CoolProp)
- Pressure: Pascals [Pa]
- Energy: Joules [J]
- Power: Watts [W]

**Display units** (in plots and output):
- Power: kilowatts [kW] = W / 1000
- Pressure: kilopascals [kPa] = Pa / 1000
- Mass flow: grams/second [g/s] = kg/s × 1000
- Time: minutes [min] = s / 60

### Important Code Patterns

#### 1. CoolProp Property Calls

Always convert temperature to Kelvin and specify quality:
```matlab
% Good
P_sat = py.CoolProp.CoolProp.PropsSI('P', 'T', T_celsius + 273.15, 'Q', 0, 'R717');

% Q = 0: saturated liquid
% Q = 1: saturated vapor
% Q = 0.5: two-phase (for saturation properties)
```

#### 2. Energy Balance Structure

```matlab
% Sequential cascade (Lines 248-402)
Q_input_net = Q_input_total * energy_factor - Q_losses
Q_NH3_circulation = NH3_flow_actual * energy_per_kg_NH3
Q_heat_pump_delivered = min(Q_heat_pump_demand, Q_NH3_circulation)
Q_water_preheat = max(0, Q_NH3_circulation - Q_heat_pump_delivered)
Q_excess = Q_input_net - Q_NH3_circulation
```

#### 3. Mass Conservation

```matlab
% Update masses (Lines 456-471)
massEvaporated = evapRate * timeStep
state.totalVaporMass = state.initialVaporMass + state.totalEvaporated
state.liquidMass = state.liquidMass - massEvaporated

% Enforce bounds
state.liquidMass = max(0, state.liquidMass)
state.totalVaporMass = max(1e-6, state.totalVaporMass)
```

#### 4. Hysteresis Control

```matlab
% Heat pump hysteresis (Lines 281-287)
if ~state.heat_pump_active && state.temperature >= params.T_heat_pump_start
    state.heat_pump_active = true;  % Turn ON
elseif state.heat_pump_active && state.temperature <= params.T_heat_pump_stop
    state.heat_pump_active = false; % Turn OFF
end
```

## Key Fixes and Features

The code includes several important fixes documented in comments:

### ✅ Fixed Issues

1. **Heat pump hysteresis** (Lines 280-287)
   - Proper start/stop thresholds
   - Prevents rapid on/off cycling

2. **Mass balance for condensation** (Lines 422-426)
   - Handles negative Q_excess (cooling)
   - Allows condensation when heat pump demand > input

3. **Water preheating priority** (Lines 340-356)
   - Heat pump gets priority
   - Preheating only gets excess after HP satisfied

4. **Improved convergence** (Lines 714-777)
   - Adaptive step size in mass balance solver
   - Better handling of edge cases

5. **Energy-driven NH3 flow logic** (Lines 306-377)
   - Flow limited by available energy
   - Not by heat pump demand

### Known Limitations

1. **Commented code**: Lines 69 and 237-245 contain incomplete Cooper correlation and old energy factor logic
2. **Simplified heat transfer**: U_overall = 10 W/m²·K is a placeholder (Line 815)
3. **No tube-side pressure drop**: Simplified fluid mechanics
4. **Fixed energy factors**: 30% (heating) and 81.3% (circulation) are empirical

## Debugging and Diagnostics

### Debug Output

The code includes extensive debugging output (controlled by `mod(step, 200) == 0`):

```matlab
% Example debug section (Lines 541-609)
if mod(step, 200) == 0
    fprintf('\n=== Step %d: Sequential Energy Cascade ===\n', step);
    fprintf('  Temperature: %.1f°C (target: %.1f°C)\n', ...);
    % ... detailed energy flow analysis
end
```

### Energy Balance Verification

```matlab
% Lines 403-411
energy_balance_error = Q_total_input - Q_total_output
if abs(energy_balance_error) > params.energy_balance_tolerance
    warning('Energy balance error: %.1f kW at step %d', ...);
end
```

Typical tolerance: 1000 W (1 kW)

### Safety Checks

```matlab
% Lines 611-622
if state.temperature >= 123
    warning('Approaching critical temperature');
    break;
end

if state.liquidMass <= 0.001
    fprintf('All liquid evaporated');
    break;
end
```

## Modifying the Model

### Common Modifications

#### 1. Change Heat Input

```matlab
% Line 18
params.heatInput = 1500;  % W/m² (increased from 1000)
```

#### 2. Adjust Heat Pump Demand

```matlab
% Line 81-82
heat_pump_demand = @(t) 50e3 + 10e3 * sin(2*pi*t/(24*3600));
% 50 kW baseline, ±10 kW daily variation
```

#### 3. Modify System Geometry

```matlab
% Lines 20-24
params.length = 2000;              % Change depth to 2000m
params.evap_length = 0.30 * params.length;  % 30% evaporator
params.cond_length = 0.20 * params.length;  % 20% condenser
params.adiab_length = 0.50 * params.length; % 50% adiabatic
```

#### 4. Change Working Fluid

```matlab
% Line 33
params.fluid = 'Water';  % Or 'CO2', 'Propane', etc.
% WARNING: Check CoolProp compatibility and temperature ranges
```

#### 5. Adjust Simulation Duration

```matlab
% Lines 34-35
params.timeStep = 60;           % 1 minute timestep
params.totalTime = 3600 * 48;   % 48 hours
```

### Adding New Outputs

To track additional variables:

1. **Preallocate array** (after Line 159):
   ```matlab
   results.your_variable = zeros(nSteps+1, 1);
   ```

2. **Store in main loop** (after Line 516):
   ```matlab
   results.your_variable(step + 1) = calculated_value;
   ```

3. **Add to plots** (in `createUnifiedThermosyphonPlots()`)

### Testing Changes

After modifications:

1. Check CoolProp initialization succeeds
2. Monitor energy balance errors
3. Verify temperature stays within valid range
4. Check for NaN or Inf values in results
5. Review plots for physical reasonableness

## Troubleshooting

### Common Issues

#### 1. CoolProp Not Found

**Error**: `'CoolProp' Python module not found`

**Solution**:
```bash
pip install CoolProp
# In MATLAB:
pyversion('/path/to/python3')  % Point to correct Python
```

#### 2. Simulation Diverges

**Symptoms**: NaN values, temperature spikes, infinite loops

**Check**:
- Time step too large (reduce `timeStep`)
- Extreme parameters (check physical reasonableness)
- Energy balance tolerance (Line 109)
- Mass balance solver convergence (Line 767)

#### 3. Plots Not Displaying

**Issue**: Figures created but not visible

**Solution**:
```matlab
set(0, 'DefaultFigureVisible', 'on');  % Before running script
```

#### 4. Slow Execution

**Causes**:
- Too many time steps (reduce `totalTime` or increase `timeStep`)
- Debug output every step (change `mod(step, 200)` to larger value)
- CoolProp calls in tight loops (consider property caching)

### Performance Optimization

1. **Reduce debug output**: Change `mod(step, 200)` to `mod(step, 500)`
2. **Increase time step**: 60s or 120s instead of 30s
3. **Disable some plots**: Comment out `create_vapor_column_plots()`
4. **Profile code**:
   ```matlab
   profile on
   Energy1
   profile viewer
   ```

## Data Interpretation

### Equilibrium Criteria (Line 431)

System reaches equilibrium when ALL conditions met:
- Temperature error < 1°C
- Excess energy < 5 kW
- No flow constraints
- Heat pump active
- Heat pump load ≥ 90%

### Status Strings

| Status | Meaning |
|--------|---------|
| `ENERGY_EQUILIBRIUM` | Steady-state operation achieved |
| `TRANSIENT` | System still evolving |
| `FLOW_LIMITED_HEATING` | Max NH3 flow reached, excess heat available |
| `INSUFFICIENT_ENERGY` | Temperature high but heat pump can't activate |

### Constraint Status

| Constraint | Meaning |
|------------|---------|
| `MAX_FLOW` | NH3 flow at maximum (0.2 kg/s) |
| `MIN_FLOW` | NH3 flow at minimum (0.01 kg/s) |
| `NO_CIRCULATION` | Heating mode, no NH3 circulation |
| `NONE` | NH3 flow unconstrained |

## AI Assistant Guidelines

### When Working with This Codebase

1. **Always preserve units**: Don't mix SI and non-SI units
2. **Maintain energy balance**: Any change must preserve energy conservation
3. **Check CoolProp validity**: Ensure temperature/pressure within fluid ranges
4. **Test incrementally**: Small changes, verify before proceeding
5. **Keep debug output**: Don't remove diagnostic prints without good reason
6. **Document changes**: Add comments explaining modifications

### Things to Avoid

❌ **Don't** remove energy balance checks
❌ **Don't** modify mass conservation logic without verification
❌ **Don't** change hysteresis thresholds without understanding impact
❌ **Don't** disable safety checks (temperature, mass limits)
❌ **Don't** assume steady-state without checking equilibrium criteria

### Recommended Improvements

If asked to enhance the model:

1. **Implement Cooper correlation** (Line 69 placeholder)
2. **Variable heat transfer coefficients** (Line 815)
3. **Pressure drop in pipes**
4. **More sophisticated geothermal model**
5. **Heat pump COP calculation**
6. **Transient wall temperature tracking**
7. **Multi-component fluid mixtures**
8. **Economic optimization**

### Code Review Checklist

When reviewing changes:

- [ ] Energy balance closes (< 1 kW error)
- [ ] Mass conservation maintained
- [ ] Units consistent throughout
- [ ] CoolProp calls valid
- [ ] Hysteresis logic correct
- [ ] Debug output informative
- [ ] Plots reflect new variables
- [ ] Results physically reasonable
- [ ] No infinite loops possible
- [ ] Safety checks in place

## Version Control

### Current State

- **Branch**: `claude/claude-md-mhxwo5vz7pwhqxqj-01Y6eDZXJzWhUrYrCvvzGN7u`
- **Latest commit**: `27fb306 Add files via upload`
- **Status**: Clean working directory

### Making Changes

When modifying the code:

```bash
# Make your changes to Energy1.m

# Commit with descriptive message
git add Energy1.m
git commit -m "Your descriptive commit message"

# Push to the feature branch
git push -u origin claude/claude-md-mhxwo5vz7pwhqxqj-01Y6eDZXJzWhUrYrCvvzGN7u
```

### Commit Message Guidelines

Good commit messages:
- `Fix heat pump hysteresis logic for proper cooling`
- `Add pressure drop calculation in evaporator section`
- `Optimize CoolProp calls to improve performance`
- `Update geothermal temperature profile with new data`

Poor commit messages:
- `update`
- `fix bug`
- `changes`

## Additional Resources

### External Dependencies

- **CoolProp**: http://www.coolprop.org/
  - Property library for fluids
  - Python interface required
  - Documentation: http://www.coolprop.org/coolprop/HighLevelAPI.html

### Thermodynamic References

- Ammonia (R717) properties: NIST Chemistry WebBook
- Heat pipe theory: Faghri's "Heat Pipe Science and Technology"
- Two-phase flow: Collier & Thome "Convective Boiling and Condensation"

### MATLAB Documentation

- Python integration: `help py`
- Plotting: `help plot`, `help subplot`
- Structures: `help struct`
- Anonymous functions: `help function_handle`

## Contact and Support

For questions about this codebase:
1. Review this CLAUDE.md file
2. Check debug output and error messages
3. Verify CoolProp installation and configuration
4. Review energy balance diagnostics

## Changelog

### [Current] - 2025-11-13
- Created comprehensive CLAUDE.md documentation
- Documented all key functions, parameters, and workflows
- Added development guidelines and troubleshooting section

### [Previous] - Unknown date
- Initial upload of Energy1.m
- Includes fixes for:
  - Heat pump hysteresis
  - Mass balance condensation
  - Water preheating priority
  - Improved convergence

---

**Note**: This is a research/development code. Use appropriate validation before applying to real engineering systems.

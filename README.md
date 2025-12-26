# Mechanistic Flow Boiling Model for Curved Downward-Facing Surfaces

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A hybrid mechanistic model for predicting heat flux during subcooled flow boiling on curved, downward-facing heated surfaces. This model is specifically designed for geometries where buoyancy acts perpendicular to the flow and toward the wall, causing bubbles to slide along the surface rather than lifting off immediately.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Physics & Model Description](#physics--model-description)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Parameters](#input-parameters)
- [Output](#output)
- [Examples](#examples)
- [Model Validation](#model-validation)
- [References](#references)
- [License](#license)

---

## Overview

Traditional flow boiling correlations (Chen, RPI Kurul-Podowski) were developed for vertical geometries where buoyancy assists bubble removal. This model addresses the unique physics of **downward-facing curved surfaces** where:

- Buoyancy acts perpendicular to flow (toward the wall)
- Bubbles slide for extended distances (10-50 mm) before lifting off
- Persistent "vapor tracks" insulate portions of the heated surface
- Surface roughness and wettability significantly affect nucleation

### Key Features

✅ **Mechanistic force balance** for bubble departure and lift-off  
✅ **Position-dependent** heat transfer along curved surface  
✅ **Surface property effects** (roughness, contact angle)  
✅ **Phase exclusion principle** for area fraction  
✅ **Easy-to-use** configuration at the top of the script  
✅ **Excel output** with pivot tables for analysis  

---

## Physics & Model Description

### Heat Flux Formulation

The total heat flux is calculated as:

```
q'' = (1 - A_f) × [F × h_fc × ΔT_sub + S × ψ × ΔT_sup^1.24 × ΔP^0.75]
```

Where:
- `A_f` = Area fraction covered by vapor (phase exclusion)
- `F` = Two-phase enhancement factor
- `h_fc` = Single-phase convective heat transfer coefficient
- `S = S_flow × S_sub` = Total suppression factor
- `ψ` = Forster-Zuber pool boiling coefficient

### Suppression Factor

The flow suppression factor is **mechanistically defined** as:

```
S_flow = D_dep / D_lift
```

Where:
- `D_dep` = Bubble departure diameter (from tangential force balance)
- `D_lift` = Bubble lift-off diameter (from radial force balance)

On upward-facing surfaces, `D_dep ≈ D_lift` → `S_flow ≈ 1`  
On downward-facing surfaces, `D_dep << D_lift` → `S_flow << 1` (significant suppression)

### Force Balance

**Tangential direction** (determines departure):
```
F_drag - F_surface_tension_t + F_buoyancy × sin(θ) - F_unsteady × sin(α) = 0
```

**Radial direction** (determines lift-off):
```
F_shear_lift + F_contact_pressure - F_surface_tension_r - F_hydro - F_buoyancy × cos(θ) - F_unsteady × cos(α) = 0
```

### Nucleation Site Density

Uses the **Benjamin-Balakrishnan correlation**, which accounts for:
- Surface roughness (Ra)
- Heater thermal properties (k, ρ, c)
- Fluid Prandtl number
- Wall superheat

```
N_a = 218.8 × Pr^1.63 × ΔT³ × γ^(-1) × Θ^(-0.4)
```

---

## Installation

### Prerequisites

- Python 3.8 or higher
- pip (Python package installer)

### Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/FlowBoilingModel.git
cd FlowBoilingModel
```

### Install Dependencies

```bash
pip install -r requirements.txt
```

Or install manually:

```bash
pip install numpy pandas pyXSteam openpyxl
```

---

## Quick Start

1. **Open** `flow_boiling_model.py` in your preferred editor

2. **Modify** the input parameters at the top of the file:

```python
# --- Operating Conditions ---
T_liquid = 95                           # Bulk liquid temperature [°C]
MASS_FLUX = 112                         # Mass flux G [kg/m²s]

# --- Superheat Range ---
SUPERHEAT_MIN = 1                       # Minimum wall superheat [°C]
SUPERHEAT_MAX = 20                      # Maximum wall superheat [°C]
SUPERHEAT_STEP = 1                      # Superheat increment [°C]

# --- Angular Positions ---
THETA_POSITIONS = [5, 10, 15]  # Positions on heater [degrees]

# --- Surface Properties ---
SURFACE_ROUGHNESS = 2.5e-6              # Surface roughness Ra [m]
CONTACT_ANGLE = 71.3                    # Static contact angle [degrees]
```

3. **Run** the model:

```bash
python flow_boiling_model.py
```

4. **Save** the results using the file dialog that appears

---

## Input Parameters

### Operating Conditions

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Liquid Temperature | `T_liquid` | °C | Bulk liquid temperature |
| Mass Flux | `MASS_FLUX` | kg/m²s | Mass flow rate per unit area |

### Superheat Range

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Minimum Superheat | `SUPERHEAT_MIN` | °C | Starting superheat value |
| Maximum Superheat | `SUPERHEAT_MAX` | °C | Ending superheat value |
| Superheat Step | `SUPERHEAT_STEP` | °C | Increment between values |

### Angular Positions

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Theta Positions | `THETA_POSITIONS` | degrees | List of angular positions on the curved heater (0° = top, 90° = bottom tangent) |

### Surface Properties

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Surface Roughness | `SURFACE_ROUGHNESS` | m | Arithmetic mean roughness Ra |
| Contact Angle | `CONTACT_ANGLE` | degrees | Static equilibrium contact angle |

### Heater Material Properties

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Thermal Conductivity | `HEATER_THERMAL_CONDUCTIVITY` | W/m·K | Heater material k_w |
| Density | `HEATER_DENSITY` | kg/m³ | Heater material ρ_w |
| Specific Heat | `HEATER_SPECIFIC_HEAT` | J/kg·K | Heater material c_w |

### Geometry

| Parameter | Variable | Units | Description |
|-----------|----------|-------|-------------|
| Heater Radius | `HEATER_RADIUS` | m | Radius of curvature |
| Heater Width | `HEATER_WIDTH` | m | Spanwise width of heater |
| Theta Start | `THETA_START` | degrees | Start of heated section |
| Theta End | `THETA_END` | degrees | End of heated section |

### Advanced Model Parameters

| Parameter | Variable | Default | Description |
|-----------|----------|---------|-------------|
| Bubble Influence Factor | `BUBBLE_INFLUENCE_FACTOR` | 0.7 | K_inf in area fraction calculation |
| Departure Time | `DEPARTURE_TIME` | 1000 ms | Characteristic bubble growth time |
| Unsteady Force Angle | `UNSTEADY_FORCE_ANGLE` | 20° | Direction of unsteady drag component |
| Lift Coefficient | `LIFT_COEFFICIENT` | 2.61 | Shear lift coefficient C_L |

---

## Output

The model generates an Excel file with the following sheets:

### Sheet 1: All Results

Complete dataset with all calculated parameters:

| Column | Units | Description |
|--------|-------|-------------|
| θ [deg] | degrees | Angular position |
| ΔT_sup [°C] | °C | Wall superheat |
| Heat Flux [kW/m²] | kW/m² | Predicted heat flux |
| D_dep [mm] | mm | Bubble departure diameter |
| D_lift [mm] | mm | Bubble lift-off diameter |
| S_flow [-] | - | Flow suppression factor |
| S_sub [-] | - | Subcooling suppression factor |
| S_total [-] | - | Total suppression factor |
| F [-] | - | Two-phase enhancement factor |
| A_f [-] | - | Vapor area fraction |
| N_a [sites/m²] | sites/m² | Nucleation site density |
| h_fc [W/m²K] | W/m²K | Convective heat transfer coefficient |

### Sheet 2: Heat Flux Summary

Pivot table with superheat as rows and angular position as columns.

### Sheet 3: S_flow Summary

Pivot table showing flow suppression factor variation.

### Sheet 4: Input Conditions

Record of all input parameters used for the calculation.

---

## Examples

### Example 1: Standard Water at Atmospheric Pressure

```python
T_liquid = 95                           # 5°C subcooling below 100°C
MASS_FLUX = 100                         # Moderate flow rate
THETA_POSITIONS = [15, 30, 45, 60]      # Four positions
SUPERHEAT_MIN = 1
SUPERHEAT_MAX = 20
SURFACE_ROUGHNESS = 1.0e-6              # 1 µm (smooth surface)
CONTACT_ANGLE = 60                      # Hydrophilic
```

### Example 2: Rough Oxidized Surface

```python
SURFACE_ROUGHNESS = 3.0e-6              # 3 µm (rough/oxidized)
CONTACT_ANGLE = 75                      # Higher contact angle due to oxide
```

### Example 3: High Flow Rate Study

```python
MASS_FLUX = 200                         # Higher mass flux
THETA_POSITIONS = np.arange(10, 80, 5).tolist()  # Fine angular resolution
```

---

## Model Validation

This model has been validated against experimental data from subcooled flow boiling experiments on curved downward-facing surfaces. Key validation points:

- Heat flux predictions within ±20% of experimental data
- Correct capture of position-dependent effects
- Accurate representation of surface roughness influence (rough surfaces show 30-40% higher heat flux)

---

## References

1. **Benjamin, R.J. and Balakrishnan, A.R.** (1996) "Nucleation site density in pool boiling of saturated pure liquids: Effect of surface microroughness and surface and liquid physical properties," *Experimental Thermal and Fluid Science*, 15(1), pp. 32-42.

2. **Mikic, B.B. and Rohsenow, W.M.** (1969) "A new correlation of pool-boiling data including the effect of heating surface characteristics," *Journal of Heat Transfer*, 91(2), pp. 245-250.

3. **Chen, J.C.** (1966) "Correlation for boiling heat transfer to saturated fluids in convective flow," *Industrial & Engineering Chemistry Process Design and Development*, 5(3), pp. 322-329.

4. **Kurul, N. and Podowski, M.Z.** (1991) "On the modeling of multidimensional effects in boiling channels," *Proceedings of the 27th National Heat Transfer Conference*, Minneapolis, MN.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## Contact

For questions or collaboration opportunities, please open an issue on GitHub.

"""
================================================================================
MECHANISTIC FLOW BOILING MODEL FOR CURVED DOWNWARD-FACING SURFACES
================================================================================

A hybrid mechanistic model for predicting heat flux during subcooled flow boiling
on curved, downward-facing heated surfaces.

This model incorporates:
- Bubble dynamics with force balance (tangential and radial)
- Benjamin-Balakrishnan nucleation site density correlation
- Surface roughness and wettability effects
- Position-dependent thermal boundary layer development
- Phase exclusion principle for area fraction
- Mechanistic suppression factor S_flow = D_dep / D_lift

Author: Sameer Osman
Institution: Khalifa University
License: MIT

================================================================================
"""

import numpy as np
import pandas as pd
from pyXSteam.XSteam import XSteam
import tkinter as tk
from tkinter import filedialog
import warnings
warnings.filterwarnings('ignore')

# ==========================================
# USER INPUT SECTION - MODIFY THESE VALUES
# ==========================================

# --- Operating Conditions ---
T_liquid = 95                           # Bulk liquid temperature [°C]
MASS_FLUX = 112                         # Mass flux G [kg/m²s]

# --- Superheat Range ---
SUPERHEAT_MIN = 1                       # Minimum wall superheat [°C]
SUPERHEAT_MAX = 20                      # Maximum wall superheat [°C]
SUPERHEAT_STEP = 1                      # Superheat increment [°C]

# --- Angular Positions ---
# Specify the theta positions (angle from top) on the curved heater [degrees]
# Example: [15, 30, 45, 60, 75] or np.arange(15, 90, 15)
THETA_POSITIONS = [5, 10, 15]  # [degrees]

# --- Surface Properties ---
SURFACE_ROUGHNESS = 2.5e-6              # Surface roughness Ra [m] (e.g., 2.5e-6 = 2.5 µm)
CONTACT_ANGLE = 71.3                    # Static contact angle [degrees]

# --- Heater Material Properties ---
HEATER_THERMAL_CONDUCTIVITY = 15.0      # Heater thermal conductivity k_w [W/m·K]
HEATER_DENSITY = 8000.0                 # Heater density ρ_w [kg/m³]
HEATER_SPECIFIC_HEAT = 500.0            # Heater specific heat c_w [J/kg·K]

# --- Geometry ---
HEATER_RADIUS = 0.523                   # Heater radius of curvature R [m]
HEATER_WIDTH = 0.02                     # Heater width (spanwise) w [m]
THETA_START = 0                         # Start of heated section [degrees]
THETA_END = 90                          # End of heated section [degrees]

# --- Model Parameters (Advanced) ---
BUBBLE_INFLUENCE_FACTOR = 0.7           # K_inf: Bubble influence area factor [-]
DEPARTURE_TIME = 1000                   # Characteristic departure time t_d [ms]
UNSTEADY_FORCE_ANGLE = 20               # Unsteady force direction α [degrees]
LIFT_COEFFICIENT = 2.61                 # Shear lift coefficient C_L [-]

# ==========================================
# PHYSICAL CONSTANTS
# ==========================================
g = 9.81                                # Gravitational acceleration [m/s²]

# ==========================================
# DO NOT MODIFY BELOW THIS LINE
# ==========================================

# Convert temperatures to Kelvin
T_l = T_liquid + 273.15

# Assign to module-level variables for functions
k_w = HEATER_THERMAL_CONDUCTIVITY
rho_w = HEATER_DENSITY
c_w = HEATER_SPECIFIC_HEAT
Ra = SURFACE_ROUGHNESS
Beta = CONTACT_ANGLE
R = HEATER_RADIUS
w = HEATER_WIDTH
theta_0 = THETA_START
theta_f = THETA_END
K_INF = BUBBLE_INFLUENCE_FACTOR
t_d = DEPARTURE_TIME
alpha = UNSTEADY_FORCE_ANGLE
CL = LIFT_COEFFICIENT

# Water properties at T_l (fitted correlations valid 60-100°C)
rho_l = 746.025 + 1.93017*T_l - 0.0036547*T_l**2                            # [kg/m³]
c_l = 9850.69 - 48.6714*T_l + 0.13736*T_l**2 - 0.000127*T_l**3              # [J/kg·K]
mu_l = 0.116947 - 0.0010053*T_l + 2.903e-6*T_l**2 - 2.806e-9*T_l**3         # [Pa·s]
k_l = -0.710696 + 0.0071857*T_l - 9.298e-6*T_l**2                           # [W/m·K]
sigma = (95130 - 0.6957*T_l - 0.2582*T_l**2) / 1000000                      # [N/m]
th_diff = k_l / rho_l / c_l                                                  # [m²/s]
Pr_l = (c_l * mu_l) / k_l                                                    # [-]

# Initialize steam table
steamTable = XSteam(XSteam.UNIT_SYSTEM_BARE)


def get_steam_properties(theta, T_w):
    """
    Get steam properties at given angular position and wall temperature.
    
    Parameters:
        theta: Angular position on heater [degrees]
        T_w: Wall temperature [K]
    
    Returns:
        P: Local pressure [kPa]
        T_sat: Saturation temperature [K]
        h_v, h_l: Vapor and liquid enthalpies [J/kg]
        hfg: Latent heat [J/kg]
        rho_g: Vapor density [kg/m³]
        mu_g: Vapor viscosity [Pa·s]
    """
    P = 101.35 + (rho_l * g * (R - R*(1 - np.cos(np.radians(theta))))) / 1000
    T_sat = steamTable.tsat_p(P/1000)
    h_v = steamTable.hV_t(T_w)
    h_l_val = steamTable.hL_t(T_w)
    hfg = h_v - h_l_val
    rho_g = steamTable.rhoV_t(T_w)
    mu_g = steamTable.my_pt(P/1000, T_w)
    return P, T_sat, h_v, h_l_val, hfg, rho_g, mu_g


def beta_func(sliding_vel, T_sat_local):
    """
    Calculate dynamic advancing and receding contact angles using 
    molecular-kinetic theory and Cox-Voinov model.
    
    Parameters:
        sliding_vel: Bubble sliding velocity [m/s]
        T_sat_local: Local saturation temperature [K]
    
    Returns:
        Beta_A: Advancing contact angle [rad]
        Beta_R: Receding contact angle [rad]
    """
    Ca = np.maximum(mu_l*sliding_vel/sigma, 1e-12)
    kb = 1.380649e-23  # Boltzmann constant
    lam = 2.5e-9       # Molecular jump length
    kappa = 3e3        # Jump frequency
    L = 1.5e-3         # Macroscopic length scale
    Ls = np.sqrt(sigma/rho_l/g)*Ca**(1/3)
    
    Beta_0_A = np.arccos((np.cos(np.radians(Beta))-(2*kb*T_sat_local/sigma/lam**2)*np.arcsinh(sliding_vel/2/kappa/lam)))
    Beta_0_R = np.arccos((np.cos(np.radians(Beta))+(2*kb*T_sat_local/sigma/lam**2)*np.arcsinh(sliding_vel/2/kappa/lam)))
    Beta_A = (Beta_0_A**3+9*Ca*np.log(L/Ls))**(1/3)
    Beta_R = (Beta_0_R**3-9*Ca*np.log(L/Ls))**(1/3)
    return Beta_A, Beta_R


def calculate_forces(r, drdt, d2rdt2, rho_g, u_eff, Beta_A, Beta_R, theta):
    """
    Calculate all forces acting on a bubble attached to the surface.
    
    Forces included:
        - Unsteady drag (bubble growth inertia)
        - Quasi-steady drag (flow resistance)
        - Shear lift (perpendicular to flow)
        - Surface tension (tangential and radial components)
        - Buoyancy
        - Contact pressure
        - Hydrodynamic pressure
    
    Parameters:
        r: Bubble radius [m]
        drdt: Bubble growth rate [m/s]
        d2rdt2: Bubble growth acceleration [m/s²]
        rho_g: Vapor density [kg/m³]
        u_eff: Effective velocity at bubble center [m/s]
        Beta_A, Beta_R: Dynamic contact angles [rad]
        theta: Angular position [degrees]
    
    Returns:
        Tuple of all force components [N]
    """
    # Unsteady Drag
    Fdu = rho_l * np.pi * r**2 * (1.5 * drdt**2 - r * d2rdt2)
    
    # Quasi-Steady Drag
    Re_b = np.abs(rho_l * u_eff * 2 * r / mu_l)
    Eo = (rho_l - rho_g) * g * (2*r)**2 / sigma
    CD1 = np.where(Re_b < 1000, 24/Re_b * (1+0.15*Re_b**0.687), 0.44)
    CD2 = np.maximum(np.minimum(16/Re_b * (1+0.15*Re_b**0.687), 48/Re_b), 8/3*Eo/(Eo+4))
    CD = (CD1 + CD2) / 2
    Fd = 0.5 * rho_l * u_eff**2 * np.pi * r**2 * CD * np.sign(u_eff) * 2
    
    # Shear Lift
    Fsl = 0.5 * rho_l * u_eff**2 * np.pi * r**2 * CL
    
    # Surface Tension
    dw = r / 7.5  # Contact diameter
    Fstt = 1.25 * dw * sigma * (np.pi * (Beta_A - Beta_R) / (np.pi**2 - (Beta_A - Beta_R)**2)) * (np.sin(Beta_A) + np.sin(Beta_R))
    Fstr = dw * sigma * (np.pi / (Beta_A - Beta_R)) * (np.cos(Beta_R) - np.cos(Beta_A))
    
    # Buoyancy
    Fb = (4/3) * np.pi * r**3 * (rho_l - rho_g) * g
    
    # Contact Pressure
    Fcp = np.pi/4 * dw**2 * 2 * sigma / (5 * r)
    
    # Hydrodynamic Pressure
    Fh = (9/8) * rho_l * u_eff**2 * (np.pi/4) * dw**2
    
    return Fdu, Fd, Fsl, Fstt, Fstr, Fb, Fcp, Fh


def get_Nsd_Benjamin(Tw_local, Tsat_local, P_local_kPa, sigma_local):
    """
    Calculate nucleation site density using Benjamin-Balakrishnan correlation.
    
    This correlation accounts for:
        - Surface roughness (Ra)
        - Heater thermal properties
        - Fluid properties
        - Wall superheat
    
    Parameters:
        Tw_local: Wall temperature [K]
        Tsat_local: Saturation temperature [K]
        P_local_kPa: Local pressure [kPa]
        sigma_local: Surface tension [N/m]
    
    Returns:
        Na: Nucleation site density [sites/m²]
        Theta_surf_param: Surface parameter [-]
    """
    gamma = np.sqrt((k_w * rho_w * c_w) / (k_l * rho_l * c_l))
    P_Pa = P_local_kPa * 1000 
    dim_group = (Ra * P_Pa) / sigma_local
    Theta_surf_param = 14.5 - 4.5 * dim_group + 0.4 * dim_group**2
    d_T = Tw_local - Tsat_local
    Na = 218.8 * (Pr_l**1.63) * (d_T**3) * (1/gamma) * (Theta_surf_param**-0.4)
    return Na, Theta_surf_param


def calculate_heat_flux(theta_deg, delta_T_superheat, mdot_flux):
    """
    Main function to calculate heat flux at given conditions.
    
    The model uses a two-step force balance approach:
    1. Tangential force balance to find departure diameter (D_dep)
    2. Radial force balance to find lift-off diameter (D_lift)
    
    Heat flux is partitioned as:
        q'' = (1 - A_f) * [F * h_fc * ΔT_sub + S * ψ * ΔT_sup^1.24 * ΔP^0.75]
    
    Parameters:
        theta_deg: Angular position on heater [degrees]
        delta_T_superheat: Wall superheat (T_w - T_sat) [°C or K]
        mdot_flux: Mass flux [kg/m²s]
    
    Returns:
        dict: Contains heat flux and all intermediate parameters
    """
    # Return structure for failed calculations
    nan_result = {
        'Heat Flux [kW/m2]': np.nan,
        'D_dep [mm]': np.nan,
        'D_lift [mm]': np.nan,
        'S_flow': np.nan,
        'S_sub': np.nan,
        'S_total': np.nan,
        'F': np.nan,
        'A_f': np.nan,
        'N_a [sites/m2]': np.nan,
        'h_fc [W/m2K]': np.nan
    }
    
    if delta_T_superheat < 0.5:
        return nan_result
    
    # Flow parameters
    u_Bulk = mdot_flux / rho_l
    Re_x = rho_l * u_Bulk * R * theta_deg * np.pi/180 / mu_l
    C_f = 0.074 * Re_x**(-1/5)
    
    # Get saturation temperature
    P_dummy, T_sat_val, _, _, _, _, _ = get_steam_properties(theta_deg, 373.15)
    T_w = T_sat_val + delta_T_superheat
    
    # Steam properties at wall temperature
    P, T_sat, h_v, h_l_val, hfg, rho_g, mu_g = get_steam_properties(theta_deg, T_w)
    
    # Bubble growth
    t = np.linspace(0.0001, 2, 10000)
    Ja = (rho_l / rho_g) * (c_l / (hfg*1000)) * (T_w - T_sat)
    C1 = np.sqrt(12/np.pi) * np.sqrt(rho_l * c_l * k_l) / (rho_g * hfg * 1000)
    Rb = (1 / 4) * C1 * (T_w - T_sat) * np.sqrt(t_d/1000)
    b_star = (2.7183 * Rb) / (C1 * (T_w - T_sat) * np.sqrt(t_d/1000))
    r = b_star * np.sqrt(12/np.pi) * Ja * np.exp(-np.sqrt(t/(t_d/1000))) * np.sqrt(th_diff * t)
    drdt = r / (2 * t) * (1 - np.sqrt(t / (t_d/1000)))
    d2rdt2 = r / (4 * t**2) * (t / (t_d/1000) - np.sqrt(t / (t_d/1000)) - 1)
    
    # Velocity profile (log-law)
    u_fr = np.maximum(u_Bulk*np.sqrt(C_f/2)*((1/0.41)*np.log(u_Bulk*np.sqrt(C_f/2)*r/mu_l*rho_l)+5.2), 0.001)
    
    # Initial sliding velocity
    u_sl = 0.00001
    u_eff = u_fr - u_sl
    Beta_A, Beta_R = beta_func(u_sl, T_sat)
    
    Fdu, Fd, Fsl, Fstt, Fstr, Fb, Fcp, Fh = calculate_forces(r, drdt, d2rdt2, rho_g, u_eff, Beta_A, Beta_R, theta_deg)
    
    # Tangential force balance -> departure diameter
    sumFt = Fd - Fstt + Fb * np.sin(theta_deg * np.pi/180) - Fdu * np.sin(alpha * np.pi/180)
    crossings = np.where(np.diff(np.sign(sumFt)) > 0)[0]
    idx_dep = crossings[0] if len(crossings) > 0 else None
    r_dep = r[idx_dep]/3 if idx_dep is not None else None
    
    if r_dep is None:
        return nan_result
    
    # Radial force balance -> lift-off diameter (iterative)
    r_lift = None
    max_iterations = 10000
    iteration = 0
    
    while r_lift is None and iteration < max_iterations:
        u_sl += 0.001
        u_eff = u_fr - u_sl
        Beta_A, Beta_R = beta_func(u_sl, T_sat)
        Fdu, Fd, Fsl, Fstt, Fstr, Fb, Fcp, Fh = calculate_forces(r, drdt, d2rdt2, rho_g, u_eff, Beta_A, Beta_R, theta_deg)
        sumFr = Fsl + Fcp - Fstr - Fh - Fb * np.cos(theta_deg * np.pi/180) - Fdu * np.cos(alpha * np.pi/180)
        crossings_fr = np.where(np.diff(np.sign(sumFr)) > 0)[0]
        idx_lift = crossings_fr[crossings_fr > idx_dep][0] if np.any(crossings_fr > idx_dep) else None
        r_lift = r[idx_lift] if idx_lift is not None else None
        iteration += 1
    
    if r_lift is None:
        return nan_result
    
    # Heat flux calculation
    Na, Theta_surf_param = get_Nsd_Benjamin(T_w, T_sat, P, sigma)
    
    D_dep = r_dep * 2
    D_lift = r_lift * 2
    A_bubble = K_INF * Na * (np.pi * D_dep**2 / 4)
    A_f = 1 - np.exp(-A_bubble)
    
    S_flow = r_dep / r_lift
    S_sub = (T_w - T_sat) / (T_w - T_l)
    S_total = S_flow * S_sub
    
    v = 0.95
    Xtt = (rho_l/rho_g)**0.5 * (mu_g/mu_l)**0.1 * v**0.9
    F = 1 + (1 / Xtt)**0.8
    
    Re_Dh = rho_l * u_Bulk * w / mu_l
    hfc = (0.2058 - 0.1396 * (theta_deg-theta_0) / (theta_f-theta_0)) * Re_Dh**0.6 * Pr_l**(1/3) * k_l / ((R*(theta_deg-theta_0)*np.pi/180*w)/(2*R*(theta_deg-theta_0)*np.pi/180 + 2*w))
    psi = 0.00122 * k_l**0.79 * c_l**0.45 * rho_l**0.49 / (sigma**0.5 * mu_l**0.29 * hfg**0.24 * rho_g**0.24)
    
    term_conv = (1-A_f) * F * hfc * (T_w - T_l)
    term_nb = (1-A_f) / (Theta_surf_param**0.4) * S_total * psi * (T_w - T_sat)**1.25 * ((steamTable.psat_t(T_w) - steamTable.psat_t(T_sat)) * 1000000)**0.75
    
    flux = term_conv + term_nb
    
    return {
        'Heat Flux [kW/m2]': flux / 1000,
        'D_dep [mm]': D_dep * 1000,
        'D_lift [mm]': D_lift * 1000,
        'S_flow': S_flow,
        'S_sub': S_sub,
        'S_total': S_total,
        'F': F,
        'A_f': A_f,
        'N_a [sites/m2]': Na,
        'h_fc [W/m2K]': hfc
    }


def main():
    """
    Main execution function.
    Runs the model for all specified conditions and saves results to Excel.
    """
    print("="*70)
    print("MECHANISTIC FLOW BOILING MODEL")
    print("For Curved Downward-Facing Surfaces")
    print("="*70)
    
    # Display input conditions
    print(f"\n--- Input Conditions ---")
    print(f"  Liquid temperature:    {T_liquid} °C")
    print(f"  Mass flux:             {MASS_FLUX} kg/m²s")
    print(f"  Superheat range:       {SUPERHEAT_MIN} - {SUPERHEAT_MAX} °C (step: {SUPERHEAT_STEP})")
    print(f"  Angular positions:     {THETA_POSITIONS} degrees")
    print(f"\n--- Surface Properties ---")
    print(f"  Surface roughness:     {SURFACE_ROUGHNESS*1e6:.2f} µm")
    print(f"  Contact angle:         {CONTACT_ANGLE}°")
    print(f"\n--- Heater Properties ---")
    print(f"  Thermal conductivity:  {HEATER_THERMAL_CONDUCTIVITY} W/m·K")
    print(f"  Radius of curvature:   {HEATER_RADIUS} m")
    
    # Generate arrays
    superheats = np.arange(SUPERHEAT_MIN, SUPERHEAT_MAX + SUPERHEAT_STEP, SUPERHEAT_STEP)
    
    # Store all results
    all_data = []
    
    print(f"\n--- Calculating ---")
    total_calcs = len(THETA_POSITIONS) * len(superheats)
    calc_count = 0
    
    for theta in THETA_POSITIONS:
        print(f"\n  θ = {theta}°:")
        for dT in superheats:
            calc_count += 1
            result = calculate_heat_flux(theta, dT, MASS_FLUX)
            
            row = {
                'θ [deg]': theta,
                'ΔT_sup [°C]': dT,
                'Heat Flux [kW/m²]': result['Heat Flux [kW/m2]'],
                'D_dep [mm]': result['D_dep [mm]'],
                'D_lift [mm]': result['D_lift [mm]'],
                'S_flow [-]': result['S_flow'],
                'S_sub [-]': result['S_sub'],
                'S_total [-]': result['S_total'],
                'F [-]': result['F'],
                'A_f [-]': result['A_f'],
                'N_a [sites/m²]': result['N_a [sites/m2]'],
                'h_fc [W/m²K]': result['h_fc [W/m2K]']
            }
            all_data.append(row)
            
            if not np.isnan(result['Heat Flux [kW/m2]']):
                print(f"    ΔT = {dT:5.1f}°C  →  q'' = {result['Heat Flux [kW/m2]']:8.2f} kW/m²")
    
    # Create DataFrame
    df = pd.DataFrame(all_data)
    
    # Count successful calculations
    valid_count = df['Heat Flux [kW/m²]'].notna().sum()
    print(f"\n--- Complete ---")
    print(f"  Successful calculations: {valid_count}/{total_calcs}")
    
    # File save dialog
    print(f"\n--- Save Results ---")
    root = tk.Tk()
    root.withdraw()
    root.attributes('-topmost', True)
    
    file_path = filedialog.asksaveasfilename(
        title="Save Results As",
        defaultextension=".xlsx",
        filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
        initialfile="flow_boiling_results"
    )
    
    if file_path:
        # Save to Excel with multiple sheets
        with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
            # Full results
            df.to_excel(writer, sheet_name='All Results', index=False)
            
            # Pivot table for heat flux
            pivot_flux = df.pivot_table(
                values='Heat Flux [kW/m²]', 
                index='ΔT_sup [°C]', 
                columns='θ [deg]'
            )
            pivot_flux.to_excel(writer, sheet_name='Heat Flux Summary')
            
            # Pivot table for S_flow
            pivot_sflow = df.pivot_table(
                values='S_flow [-]', 
                index='ΔT_sup [°C]', 
                columns='θ [deg]'
            )
            pivot_sflow.to_excel(writer, sheet_name='S_flow Summary')
            
            # Input conditions summary
            conditions = pd.DataFrame({
                'Parameter': [
                    'Liquid Temperature [°C]',
                    'Mass Flux [kg/m²s]',
                    'Surface Roughness [µm]',
                    'Contact Angle [°]',
                    'Heater Thermal Conductivity [W/m·K]',
                    'Heater Radius [m]',
                    'Heater Width [m]',
                    'Bubble Influence Factor K_inf',
                    'Lift Coefficient C_L'
                ],
                'Value': [
                    T_liquid,
                    MASS_FLUX,
                    SURFACE_ROUGHNESS * 1e6,
                    CONTACT_ANGLE,
                    HEATER_THERMAL_CONDUCTIVITY,
                    HEATER_RADIUS,
                    HEATER_WIDTH,
                    BUBBLE_INFLUENCE_FACTOR,
                    LIFT_COEFFICIENT
                ]
            })
            conditions.to_excel(writer, sheet_name='Input Conditions', index=False)
        
        print(f"  Results saved to: {file_path}")
    else:
        print("  Save cancelled.")
    
    root.destroy()
    
    print(f"\n{'='*70}")
    print("Model execution complete!")
    print(f"{'='*70}")
    
    return df


if __name__ == "__main__":
    results = main()

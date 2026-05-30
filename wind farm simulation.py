"""
Wind farm simulation for a 90 MW crescent-arc layout in Al-Tafilah, Jordan.
Implements Jensen (Park) and Gaussian (Bastankhah & Porté-Agel, 2014) wake models.
"""

import math
import numpy as np
from scipy.special import gamma
from scipy.stats import weibull_min

# ------------------------------
# 1. Parameters
# ------------------------------
D = 112.0          # rotor diameter (m)
R = D / 2.0        # rotor radius (m)
Ct = 0.80          # thrust coefficient (constant approximation)
k_jensen = 0.075   # wake decay constant for Jensen model
TI = 0.06          # ambient turbulence intensity (6% for onshore)

# Turbine coordinates (from Table 5 in report)
turbines = np.array([
    [-647.05, 2414.81], [-178.35, 2493.63], [296.80, 2482.32], [761.22, 2381.29],
    [1198.12, 2194.20], [1591.73, 1927.80], [1927.80, 1591.73], [2194.20, 1198.12],
    [2381.29, 761.22], [2482.32, 296.80], [2493.63, -178.35], [2414.81, -647.05],
    [-248.70, 1697.88], [149.56, 1709.47], [539.76, 1628.90], [900.86, 1460.52],
    [1213.40, 1213.40], [1460.52, 900.86], [1628.90, 539.76], [1709.47, 149.56],
    [1697.88, -248.70], [1594.76, -633.56], [-241.22, 900.24], [34.85, 931.35],
    [307.82, 879.70], [553.44, 749.89], [749.89, 553.44], [879.70, 307.82],
    [931.35, 34.85], [900.24, -241.22]
])

# ------------------------------
# 2. Wind resource generation (Weibull with monthly means)
# ------------------------------
monthly_speeds_45m = {
    1: 10.1, 2: 11.1, 3: 10.5, 4: 9.1, 5: 7.8, 6: 8.0,
    7: 7.8, 8: 7.6, 9: 6.0, 10: 5.9, 11: 8.9, 12: 9.0
}
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
scale_height = (85 / 45) ** 0.143   # power law exponent 0.143

np.random.seed(42)   # reproducible results

wind_speeds = []
wind_dirs = []

for month, days in enumerate(days_per_month, start=1):
    v45 = monthly_speeds_45m[month]
    v_hub = v45 * scale_height
    # Weibull scale factor to match monthly mean (shape k=2.2)
    c_month = v_hub / gamma(1 + 1/2.2)
    for _ in range(days):
        for _ in range(24):
            speed = max(0.5, np.random.weibull(2.2) * c_month)
            wind_speeds.append(speed)
            wd = (315 + np.random.normal(0, 45)) % 360
            wind_dirs.append(wd)

wind_speeds = np.array(wind_speeds)
wind_dirs = np.array(wind_dirs)

# Verify Weibull fit
k_fit, loc, c_fit = weibull_min.fit(wind_speeds, floc=0)
print(f"Annual mean wind speed: {wind_speeds.mean():.2f} m/s (target 9.27 m/s)")
print(f"Weibull fit: c={c_fit:.2f} m/s (target 10.47), k={k_fit:.2f} (target 2.20)")

# ------------------------------
# 3. Helper functions
# ------------------------------
def rotate_coords(coords, wind_dir_deg):
    """Rotate turbine coordinates so wind blows from left to right (x direction)."""
    theta = math.radians(270 - wind_dir_deg)
    c, s = math.cos(theta), math.sin(theta)
    rot = np.zeros_like(coords)
    rot[:, 0] = c * coords[:, 0] + s * coords[:, 1]
    rot[:, 1] = -s * coords[:, 0] + c * coords[:, 1]
    return rot

def power_curve(ws):
    """Vestas V112 3.0 MW power curve (simplified cubic)."""
    if ws < 3.0 or ws >= 25.0:
        return 0.0
    if ws >= 12.0:
        return 3.0
    return 3.0 * ((ws - 3.0) / (12.0 - 3.0)) ** 3

# ------------------------------
# 4. Jensen (Park) wake model
# ------------------------------
def jensen_wake_speeds(turbines_xy, u_inf, wind_dir, k=k_jensen):
    rot = rotate_coords(turbines_xy, wind_dir)
    N = len(turbines_xy)
    speeds = np.full(N, u_inf)
    for i in range(N):
        max_def = 0.0
        xi = rot[i, 0]
        yi = rot[i, 1]
        for j in range(N):
            if rot[j, 0] >= xi - 1.0:
                continue   # j is not upstream
            dx = xi - rot[j, 0]
            dy = abs(yi - rot[j, 1])
            wake_radius = R + k * dx
            if wake_radius <= 0 or dy > wake_radius:
                continue
            # Overlap area factor (partial wake)
            if dy + R <= wake_radius:
                overlap = 1.0
            else:
                # Circular segment intersection
                d1 = max(-1.0, min(1.0, (dy**2 + wake_radius**2 - R**2) / (2 * dy * wake_radius)))
                d2 = max(-1.0, min(1.0, (dy**2 + R**2 - wake_radius**2) / (2 * dy * R)))
                A1 = wake_radius**2 * math.acos(d1) - wake_radius**2 * d1 * math.sqrt(1 - d1**2)
                A2 = R**2 * math.acos(d2) - R**2 * d2 * math.sqrt(1 - d2**2)
                overlap = max(0.0, min(1.0, (A1 + A2) / (math.pi * R**2)))
            deficit = overlap * (1 - math.sqrt(1 - Ct)) * (D / (D + 2 * k * dx)) ** 2
            max_def = max(max_def, deficit)
        speeds[i] = u_inf * (1 - max_def)
    return speeds

# ------------------------------
# 5. Gaussian wake model (Bastankhah & Porté-Agel, 2014)
# ------------------------------
def gaussian_wake_speeds(turbines_xy, u_inf, wind_dir, ti=TI):
    rot = rotate_coords(turbines_xy, wind_dir)
    N = len(turbines_xy)
    speeds = np.full(N, u_inf)
    # Wake growth rate (Niayifar & Porté-Agel, 2016)
    k_star = 0.3837 * ti + 0.003678
    epsilon = 0.2 * math.sqrt(0.5 * (1 + math.sqrt(1 - Ct)) / math.sqrt(max(1e-6, 1 - Ct)))
    for i in range(N):
        sum_deficit_sq = 0.0
        xi = rot[i, 0]
        yi = rot[i, 1]
        for j in range(N):
            if rot[j, 0] >= xi - 1.0:
                continue
            dx = xi - rot[j, 0]
            dy = abs(yi - rot[j, 1])
            xD = dx / D
            if xD <= 0:
                continue
            sigma = D * (k_star * xD + epsilon)
            # Centreline velocity deficit
            C0 = 1 - math.sqrt(max(0.0, 1 - Ct / (8 * (sigma / D)**2)))
            deficit = C0 * math.exp(-0.5 * (dy / sigma)**2)
            sum_deficit_sq += deficit**2
        speeds[i] = u_inf * (1 - math.sqrt(sum_deficit_sq))
    return speeds

# ------------------------------
# 6. Annual simulation
# ------------------------------
N_hours = 8760
N_turb = len(turbines)

aep = {"No Wake": 0.0, "Jensen (Park)": 0.0, "Gaussian": 0.0}

for h in range(N_hours):
    ws = wind_speeds[h]
    wd = wind_dirs[h]
    # No wake: every turbine sees freestream
    aep["No Wake"] += N_turb * power_curve(ws)
    # Jensen
    speeds_j = jensen_wake_speeds(turbines, ws, wd)
    aep["Jensen (Park)"] += np.sum([power_curve(s) for s in speeds_j])
    # Gaussian
    speeds_g = gaussian_wake_speeds(turbines, ws, wd)
    aep["Gaussian"] += np.sum([power_curve(s) for s in speeds_g])

# Convert Wh to GWh (1 GWh = 1e9 Wh)
for key in aep:
    aep[key] /= 1e9

print("\n--- Annual Simulation Results ---")
print(f"No Wake:        {aep['No Wake']:.2f} GWh/yr")
print(f"Jensen (Park):  {aep['Jensen (Park)']:.2f} GWh/yr  |  Wake loss: {(aep['No Wake'] - aep['Jensen (Park)'])/aep['No Wake']*100:.2f}%")
print(f"Gaussian:       {aep['Gaussian']:.2f} GWh/yr  |  Wake loss: {(aep['No Wake'] - aep['Gaussian'])/aep['No Wake']*100:.2f}%")

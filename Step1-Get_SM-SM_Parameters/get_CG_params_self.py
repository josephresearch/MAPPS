import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy import integrate
from scipy.interpolate import UnivariateSpline

######### INPUT PARAMETERS #########
# INPUT SIGMA PARAMETER FOR the NEW BEAD 

def get_profile(file):
    with open(file, 'r') as f:
        data = [line.split() for line in f if not line.startswith(('#', '@'))]
        data = list(zip(*data))  # Transpose the data for plotting
        x, y = list(map(float, data[0])), list(map(float, data[1]))
        # shift the profile at the bulk region to 0  (bulk = 1.4 nm)
        closest_x = min(x, key=lambda v: abs(v - 1.4))    
        index_closest_x = x.index(closest_x)
        shift_value = y[index_closest_x]
        y = [(val - shift_value) for val in y]  # maybe needed conversion to kB*T !, kcal/mol for now 
        return x, y 

def extract_first_well(x, y):
    y = np.array(y)
    well_indices = np.where(y < 0)[0]
    start_index = well_indices[0]
    end_index = start_index + np.argmax(np.diff(well_indices) != 1)
    if end_index == start_index:  # this is for the case when there is one well (CG lammps output)
        end_index = well_indices[-1]
    return x[start_index:end_index + 1], y[start_index:end_index + 1]

def integrate_well(x_well, y_well):
    return integrate.trapezoid(y_well, x_well)

def separate_decreasing_region(x, y):
    min_index = y.index(min(y))
    return x[:min_index + 1], y[:min_index + 1], x[min_index + 1:], y[min_index + 1:]

def interpolate(xbeforemin, ybeforemin):
    spline = UnivariateSpline(xbeforemin, ybeforemin, k=3, s=0)
    x_interpolated = np.linspace(min(xbeforemin), max(xbeforemin), 500)
    y_interpolated = spline(x_interpolated)
    return x_interpolated, y_interpolated

def get_integral(x, y):
    x_well, y_well = extract_first_well(x, y)
    area_under_the_well = abs(integrate_well(x_well, y_well))
    return x_well, y_well, area_under_the_well

# compute the interaction strength for YY pair (from my simulations with AMBER03)
YY_AA_profile_x, YY_AA_profile_y = get_profile('compounds/YY/profile.xvg')
YY_int_AA = get_integral(YY_AA_profile_x, YY_AA_profile_y)[2]
epsilon_YY_AA = np.min(YY_AA_profile_y)

def find_sigma(x, y):
    min_difference = float('inf')
    sigma = None
    for x_val, y_val in zip(x, y):
        difference = abs(y_val - 0)
        if difference < min_difference:
            min_difference = difference
            sigma = x_val
    return sigma

# WF potential energy function 
def phi_WF(r_dist, eps, sigma, nu, mu):
    rc = sigma * 3.0
    alpha = 2 * nu * ((rc / sigma) ** (2 * mu)) * ((1 + 2 * nu) / (2 * nu * ((rc / sigma) ** (2 * mu) - 1))) ** (2 * nu + 1)
    phi = eps * alpha * ((sigma / r_dist) ** (2 * mu) - 1) * ((abs(rc / r_dist) ** (2 * mu) - 1)) ** (2 * nu)
    return phi

# get the area under the curve from the CG potential 
r_CG = np.linspace(4.5, 20, 1000)

def CG_well(x, y):
    y_attractive = [el for el in y if el < 0]
    x_attractive = [x[ind] for ind, el in enumerate(y) if el < 0]
    return x_attractive, y_attractive

# get the value for TYR-TYR beads interaction 
phi_WF_YY = phi_WF(r_CG, 0.419186, 6.733634028253079, 1, 2)
epsilon_YY_CG = abs(np.min(phi_WF_YY))
x_well, y_well = CG_well(r_CG, phi_WF_YY)
YY_int_CG = abs(integrate_well(x_well, y_well))

# iterate through different epsilon for a CG bead, sigma is fixed 
r = np.linspace(2.5, 20, 1000)

def find_bead_params(sigma_bead, bead_int_AA, epsilon_bead_AA_relative):
    eps_values = np.linspace(0.005, 1.5, 10000)
    best_value = None
    best_score = float('inf')  # initialize with a high value
    nu = 1 
    mu_values = [1, 2, 3, 4, 5, 6, 7]
    scores = []
    values = []
    for mu_value in mu_values:
        for value in eps_values:
            phi_bead = phi_WF(r, value, sigma_bead, nu, mu_value)
            epsilon_bead_CG = value
            val_x_well, val_y_well = CG_well(r, phi_bead)
            val_integral = abs(integrate_well(val_x_well, val_y_well)) / YY_int_CG
            epsilon_bead_CG_relative = abs(epsilon_bead_CG / epsilon_YY_CG)
            score1 = abs(epsilon_bead_CG_relative - epsilon_bead_AA_relative) / epsilon_bead_AA_relative
            score2 = abs(val_integral - bead_int_AA) / bead_int_AA
            if max(score1, score2) < best_score:
                best_score = max(score1, score2)
                best_value = value
        scores.append(best_score)
        values.append(best_value)
    best_score_of_all_mu = np.min(scores)
    result_index = scores.index(best_score_of_all_mu)
    if best_score_of_all_mu < 0.2:       
        return values[result_index], nu, mu_values[result_index], best_score_of_all_mu
    else:   
        print("No values found, try changing mu and nu parameters further")
        return [0] * 4
    
def find_parameters(files_location, new_bead_self_sigma):
    new_bead_sigma = new_bead_self_sigma
    xvg_file = files_location

    # compute the area under the attractive well 
    x, y = get_profile(xvg_file)
    before_well_x, before_well_y, remaining_part_x, remaining_part_y = separate_decreasing_region(x, y)
    x_interpolated, y_interpolated = interpolate(before_well_x, before_well_y)
    x_full = np.append(x_interpolated, remaining_part_x)
    y_full = np.append(y_interpolated, remaining_part_y) 
    x_well, y_well, area_under_the_well = get_integral(x_full, y_full) 
    bead_int_AA = area_under_the_well / YY_int_AA  # relative to YY
    epsilon_bead_AA = np.min(y_full)
    epsilon_bead_AA_relative = abs(epsilon_bead_AA / epsilon_YY_AA)

    # find epsilon for a new bead
    print("PMF:", xvg_file)
    print("new bead sigma:", "{:.3f}".format(new_bead_sigma), "[A]")
    bead_int_AA_relative = bead_int_AA  # int strength from AA profiles, relative to AA YY integral value 
    new_bead_eps, nu, mu, best_score = find_bead_params(new_bead_sigma, bead_int_AA_relative, epsilon_bead_AA_relative)
    print("new bead eps:", "{:.3f}".format(new_bead_eps), "[kcal/mol]")
    print("nu:", nu, "mu:", mu)
    print("maximum deviation:", best_score)
    print("_" * 30)
    
    return new_bead_eps, new_bead_sigma, nu, mu

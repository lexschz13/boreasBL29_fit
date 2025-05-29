import numpy as np
from scipy.special import voigt_profile as voigt, gamma as gamma_fun



def lorentz(x, x0, A, g):
    return A * (g / ((x-x0)**2 + g**2)) / np.pi #Lorentz distribution


def gauss(x, x0, A, sigma):
    return A/(np.sqrt(2*np.pi)*sigma) * np.exp(-0.5*((x-x0)/sigma)**2)


def general_gauss(x, x0, A, alpha, beta):
    return A * beta/(2*alpha*gamma_fun(1/beta)) * np.exp((np.abs(x-x0)/alpha)**beta)


def modvoigt(x, x0, A, sigma, g):
    return A*voigt(x-x0, sigma, g)


def _include_constrained_param(const_params, params_per_peak, ratio):
    i = params_per_peak+1
    return tuple(const_params[:i]) + (ratio*const_params[1],) + tuple(const_params[i:])


# def mixed_peak(x, x0, A, g, sigma, r=0.3):
#     if r<0 or r>1:
#         raise ValueError("The reason between lorentzian and gaussian must be between 0 and 1")
#     return r*lorentz(x, x0, A, g) + (1-r)*gauss(x, x0, A, sigma)
    


def multi_lorentz(x, ratio=5, nlor=3, *fitcoef):
    param = fitcoef[:4] + (ratio*fitcoef[1],) + fitcoef[4:]
    profile = np.zeros_like(x)
    for n in range(nlor):
        profile += lorentz(x, param[3*n+0], param[3*n+1], param[3*n+2])
    return profile + param[-1]


def multi_gauss(x, ratio=5, nlor=3, *fitcoef):
    param = fitcoef[:4] + (ratio*fitcoef[1],) + fitcoef[4:]
    profile = np.zeros_like(x)
    for n in range(nlor):
        profile += gauss(x, param[3*n+0], param[3*n+1], param[3*n+2])
    return profile + param[-1]


def multi_general_gauss(x, ratio=5, npeaks=3, *fitcoef):
    param = _include_constrained_param(fitcoef, 4, ratio)
    profile = np.zeros_like(x)
    for n in range(npeaks):
        profile += general_gauss(x, param[4*n+0], param[4*n+1], param[4*n+2], param[4*n+3])
    return profile + param[-1]


def multi_voigt(x, ratio=5, npeaks=3, *fitcoef):
    param = _include_constrained_param(fitcoef, 4, ratio)
    profile = np.zeros_like(x)
    for n in range(npeaks):
        profile += modvoigt(x, param[4*n+0], param[4*n+1], param[4*n+2], param[4*n+3])
    return profile + param[-1]



_params_per_peak = {"voigt":4,
                   "lorentz":3,
                   "gauss":3,
                   "gen_gauss":4}
_single_func = {"voigt":modvoigt,
               "lorentz":lorentz,
               "gauss":gauss,
               "gen_gauss":general_gauss}
_multi_func = {"voigt":multi_voigt,
              "lorentz":multi_lorentz,
              "gauss":multi_gauss,
              "gen_gauss":multi_general_gauss}

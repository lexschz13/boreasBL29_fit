import numpy as np



# These indices must be modified in function of the position of metadata on corresponding .dat lines
data_cook_position = {"energy": 1,
                      "spectra_num": 6,
                      "spectra_denom": 5,
                      "scan_name": 1,
                      "temperature": 6,
                      "measured temperature": 5,
                      "pol angle": 1,
                      "magnet": 5,
                      "incidence angle": 1,
                      "zsample": 4}


# Allowed fitting modes
_allowed_modes = ["voigt", "lorentz", "gauss", "gen_gauss"]


# Guess for fittings
guess = {"voigt": (529.118, 0.6, 0.326, 0.123,
                   529.941,     0.521, 0.203,
                   533.392, 0.8, 1.946, 0.321,
                   0),
         "gen_gauss": (529.118, 0.6, 0.326, 1.236,
                        529.941,     0.521, 1.897,
                        533.392, 0.8, 1.946, 1.456,
                        0),
         "lorentz": (529.118, 0.6, 0.326,
                     529.941,     0.521,
                     533.392, 0.8, 1.946,
                     0),
         "gauss": (529.118, 0.6, 0.326,
                   529.941,     0.521,
                   533.392, 0.8, 1.946,
                   0)}

# Bounds for fittings
limits = {"voigt": (np.array([529, 0, 0, 0,
                              529.5,  0, 0,
                              533, 0, 0, 0,
                              0]),
                    np.array([530.5, np.inf, np.inf, np.inf,
                              531.5,           np.inf, np.inf,
                              700.0, np.inf, np.inf, np.inf,
                              np.inf])),
          "gen_gauss": (np.array([529, 0, 0, 0,
                                    529.5,  0, 0,
                                    533, 0, 0, 0,
                                    0]),
                        np.array([530.5, np.inf, np.inf, np.inf,
                                    531.5,           np.inf, np.inf,
                                    700.0, np.inf, np.inf, np.inf,
                                    np.inf])),
         "lorentz": (np.array([529, 0, 0,
                               529.5,  0,
                               533, 0, 0,
                               0]),
                     np.array([530.5, np.inf, np.inf,
                               531.5,           np.inf,
                               700.0, np.inf, np.inf,
                               np.inf])),
         "gauss": (np.array([529, 0, 0,
                             529.5,  0,
                             533, 0, 0,
                             0]),
                   np.array([530.5, np.inf, np.inf,
                             531.5,           np.inf,
                             700.0, np.inf, np.inf,
                             np.inf]))}

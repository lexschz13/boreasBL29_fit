import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import csv
import os
from .profiles import (modvoigt, multi_voigt,
                        lorentz, multi_lorentz, 
                        gauss, multi_gauss,
                        general_gauss, multi_general_gauss,
                        _include_constrained_param, _params_per_peak, _single_func, _multi_func)
from .defaults import guess, limits, _allowed_modes



def _fit_peaks(fit_func, datax, datay, init_guess, bounds, ratio=5, npeaks=3):
    """
    Fit spectroscopy data with a n-peak model.

    Parameters
    ----------
    fit_func : func
        Model to fit.
    datax : numpy.ndarray,list
        Array with x-coordinates of data to fit.
    datay : numpy.ndarray,list
        Array with y-coordinates of data to fit.
    init_guess : tuple,list
        Initial guess of fitting parameters.
    bounds : tuple,list(numpy.array)
        Bounds for fitting parameters.
    ratio : float, optional
        Ratio between two first peaks. The default is 5.
    npeaks : int, optional
        Number of peaks. The default is 3.

    Returns
    -------
    popt : numpy.array
        Fitting parameters.
    pcov : numpy.array
        Covariance matrix.

    """
    # assert len(init_guess) == nlor*3 # 3 args per peak
    assert len(datax) == len(datay)
    
    def single_arg_func(x, *fitcoef): # Curve fit works assuming that first argument is variable and other ones are fitting parameters
        return fit_func(x, ratio, npeaks, *fitcoef)
    
    return curve_fit(single_arg_func, datax, datay, p0=init_guess, bounds=bounds, maxfev=10000)



def _data_from_csvline(line, emin=520, emax=560):
    """
    Takes a line from .csv and returns the energy and spectra.

    Parameters
    ----------
    line : list(string)
        Strings of a line of a .csv file.
    emin : float, optional
        Bottom of energy window. The default is 520.
    emax : float, optional
        Top of energy window. The default is 560.

    Returns
    -------
    energy : numpy.ndarray
        Energy data inside the window.
    spectra : numpy.ndarray
        Spectra data inside the window.

    """
    data = np.array(line[8:]).astype(np.float64)
    energy = data[:data.size//2]
    spectra = data[data.size//2:]
    wh = np.where((energy >= emin) * (energy <= emax))
    return energy[wh], spectra[wh]



def single_fit(csvname, scan, mode="voigt", emin=520, emax=560, ratio=5, npeaks=3, plot_peaks=True):
    """
    Makes a fitting of a single scan from name.

    Parameters
    ----------
    csvname : string
        Name of .csv file where data is stored.
    scan : string
        Name of scan to fit.
    mode : string, optional
        Kind of profile for peaks. The default is "voigt".
    emin : float, optional
        Bottom of energy window. The default is 520.
    emax : float, optional
        Top of energy window. The default is 560.
    ratio : float, optional
        Ratio between two first peaks. The default is 5.
    npeaks : int, optional
        Number of peaks. The default is 3.
    plot_peaks : bool, optional
        If True plot the fitting peaks. The default is True.

    Returns
    -------
    energy : numpy.ndarray
        Energy data inside the window.
    fitted : numpy.ndarray
        Fitted spectra.
    spectra : numpy.ndarray
        Spectra data inside the window.
    gap : float
        Gap between two first peaks.
    errgap : float
        Gap deviation.
    rel_errgap :  float
        Relative gap error.
    cost : float
        Cost of least-squares.
    rel_cost : float
        Relative cost of least-squares.
    popt : numpy.array
        Fitting parameters.
    pcov : numpy.array
        Covariance matrix.
    metadata : dictionary
        Metadata of the scan.

    """
    if not mode in _allowed_modes:
        raise ValueError("Invalid mode")
    params_per_peak = _params_per_peak[mode]
    single_func = _single_func[mode]
    multi_func = _multi_func[mode]
    pinit = guess[mode]
    bounds = limits[mode]
    flname = csvname + ".csv"
    if not flname in os.listdir():
        raise ValueError("File not found")
        return
    
    with open(flname, mode='r', newline='') as csvfl:
        reader = csv.reader(csvfl)
        for line in reader:
            if line[0]==scan:
                print("Temperature: %i\n" % int(line[1]),
                      "Polarization: %s\n" % line[3],
                      "Magnet: %s\n" % line[4],
                      "Incidence angle: %.3f\n" % float(line[5]),
                      "Sample: %s" % line[7])
                metadata = {"Temperature": int(line[1]),
                            "Polarization": line[3],
                            "Magnet": line[4],
                            "Incidence angle": float(line[5]),
                            "Sample": line[7]}
                energy, spectra = _data_from_csvline(line, emin, emax)
                if energy.size == 0:
                    raise ValueError("Window outside energy data")
                popt, pcov = _fit_peaks(multi_func, energy, spectra, pinit, bounds, npeaks=npeaks)
                # error = np.sqrt(np.diag(pcov))
                pvar = np.diag(pcov)
                fitted = multi_func(energy, ratio, npeaks, *popt)
                param = _include_constrained_param(popt, params_per_peak, ratio)
                offset = popt[-1]
                if plot_peaks:
                    fig,ax = plt.subplots()
                    ax.plot(energy, spectra-offset, lw=0.8, label="Data")
                    ax.plot(energy, fitted-offset, lw=0.8, label="Fitted")
                    for n in range(npeaks):
                        ax.plot(energy, single_func(energy, *param[params_per_peak*n:params_per_peak*(n+1)]), lw=0.8, label="Peak %i" % (n+1))
                    ax.set_xlabel("Energy (eV)")
                    ax.set_ylabel("Spectra")
                    ax.legend(loc=0)
                    ax.set_title(scan)
                    
                    gig,bx = plt.subplots()
                    bx.plot(energy, spectra-fitted, lw=0.8, color='r')
                    bx.set_xlabel("Energy (eV)")
                    bx.set_ylabel("Error")
                    bx.set_title(scan)
                
                gap = popt[params_per_peak] - popt[0]
                errgap = np.sqrt(pvar[params_per_peak]+pvar[0] - 2*pcov[0,params_per_peak])
                cost = np.sqrt(np.sum((fitted-spectra)**2))
                return (energy, # Energy points
                        fitted, # Fitted data
                        spectra, # Original data
                        gap, # Gap
                        errgap, # Gap error
                        errgap / gap, # Relative error gap
                        cost, # Cost
                        cost / np.sqrt(np.sum(spectra**2)), # Relative cost
                        popt, # Fitting parameters
                        pcov, # Covariance
                        metadata, # Metadata of the scan
                        ) 
    
    raise ValueError("Non-valid scan")
    return


def compute_gap(csvname, pol, sample, B, incidence_angle,
                mode="voigt", emin=520, emax=560, ratio=5, npeaks=3,
                plot_gaps=True, ylim_gaps=(0,1)):
    """
    For a sample, a polarization, a switch of magnetic field and an incidence angle fits the gap in funciton of temperture.

    Parameters
    ----------
    csvname : string
        Name of .csv file where data is stored.
    pol : string
        Polarization.
    sample : string
        Name of the sample.
    B : string
        Magnetic field "On" of "Off".
    incidence_angle : int
        Incidence angle of the beam.
    mode : string, optional
        Kind of profile for peaks. The default is "voigt".
    emin : float, optional
        Bottom of energy window. The default is 520.
    emax : float, optional
        Top of energy window. The default is 560.
    ratio : float, optional
        Ratio between two first peaks. The default is 5.
    npeaks : int, optional
        Number of peaks. The default is 3.
    plot_gaps : bool, optional
        If True plot the gaps in function of temperature. The default is True.
    ylim_gaps : tuple(float), optional
        Plotting window of gaps. The default is (0,1).

    Returns
    -------
    T : numpy.ndarray
        Temperatures.
    gap : numpy.ndarray
        Gap.
    errorbar : numpy.ndarray
        Errors of the gaps.
    rel_errgap :  numpy.ndarray
        Relative gap error.
    cost : numpy.ndarray
        Cost of least-squares.
    rel_cost : numpy.ndarray
        Relative cost of least-squares.

    """
    if not mode in _allowed_modes:
        raise ValueError("Invalid mode")
    params_per_peak = _params_per_peak[mode]
    # single_func = _single_func[mode]
    multi_func = _multi_func[mode]
    pinit = guess[mode]
    bounds = limits[mode]
    flname = csvname + ".csv"
    if not flname in os.listdir():
        raise ValueError("File not found")
        return
    
    avg_spectra = {}
    counts = {}
    with open(flname, mode='r', newline='') as csvfl:
        reader = csv.reader(csvfl)
        for l,line in enumerate(reader):
            if l!=0:
                if pol==line[3] and B==line[4] and incidence_angle==int(float(line[5])) and sample==line[7]:
                    T = int(line[1])
                    energy, spectra = _data_from_csvline(line, emin, emax)
                    if energy.size == 0:
                        continue
                    print(line[0], "%iK" % T)
                    if T in avg_spectra.keys():
                        avg_spectra[T] = (avg_spectra[T]*counts[T] + spectra) / (counts[T] + 1)
                        counts[T] += 1
                    else:
                        avg_spectra[T] = spectra
                        counts[T] = 1
    
    Tls = np.sort(list(avg_spectra.keys()))
    gaps = []
    error_bars = []
    cost = []
    rel_cost = []
    for T in Tls:
        popt, pcov = _fit_peaks(multi_func, energy, avg_spectra[T], pinit, bounds, npeaks=npeaks)
        gaps.append(popt[params_per_peak] - popt[0])
        error_bars.append(np.sqrt(pcov[0,0] + pcov[params_per_peak,params_per_peak] - 2*pcov[0,params_per_peak]))
        fitted = multi_func(energy, ratio, npeaks, *popt)
        cost.append(np.sqrt(np.sum((avg_spectra[T]-fitted)**2)))
        rel_cost.append(np.sqrt(np.sum((avg_spectra[T]-fitted)**2)) / np.sqrt(np.sum(avg_spectra[T]**2)))
    
    Tarr = np.array(Tls)
    gaps = np.array(gaps)
    error_bars = np.array(error_bars)
    cost = np.array(cost)
    rel_cost = np.array(rel_cost)
    
    if plot_gaps:
        fig, ax = plt.subplots()
        ax.errorbar(Tarr, gaps, yerr=error_bars, lw=0.8, marker='o', ms=5, ecolor='k', color='b')
        ax.set_ylim(*ylim_gaps)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Gap (eV)")
        ax.set_title("%s %s magnet %s incidence %i %s fit" % (sample,
                                                            pol,
                                                            B,
                                                            incidence_angle,
                                                            mode))
        
        gig, bx = plt.subplots()
        bx.plot(Tarr, rel_cost*100, lw=0.8, marker='o', ms=5, color='r')
        bx.set_xlabel("Temperature (K)")
        bx.set_ylabel("Cost (%)")
        bx.set_title("%s %s magnet %s incidence %i %s fit" % (sample,
                                                            pol,
                                                            B,
                                                            incidence_angle,
                                                            mode))
    
    
    return (Tarr, # Temperature
            gaps, # Gaps
            error_bars, # Error of gaps
            error_bars / gaps, # Relative gaps errors
            cost, # Cost
            rel_cost, # Relative cost
            )
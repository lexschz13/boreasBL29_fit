import numpy as np
import csv
import os
from .defaults import data_cook_position as dcp





def _read_boreas_file(filename, disc_mode="region", sample_list=[], region_limits=[], exceptions=([],)):
    """
    Reads a .dat file from boreas BL29 and convert the usefull metadata and cooked data into a datasets.

    Parameters
    ----------
    filename : string
        Name of the .dat file from experimental data and metadata is extracted.
    disc_mode : string, optional
        Strategy for sample distinguishment from position. The default is "region".
    sample_list : list(string), optional
        Names of the samples. If list is empty there will not be any discrimination and all samples will recive the name "no-name". The default is [].
    region_limits : list(int,float), optional
        Limits of the regions occupied by samples with temperature changes. There must be one less limit than samples. The default is [].
    exceptions : tuple(list(int))
        If some sample goes out of range you can put an exception. The default is ([],).

    Returns
    -------
    data_sets : dictionary
        This storages all usefull metadata, energy and spectra.

    """
    data_sets = []
    append_new = True
    fl = open(filename + ".dat")
    lines = [i.strip().split() for i in fl.readlines()]
    data = []
    for line in lines:
        if append_new:
            data_arr = np.array(data).astype(np.float64)
            if data_arr.ndim >= 2:
                # Cooking data
                data_sets[-1]["Energy"] = data_arr[:,dcp["energy"]]
                data_sets[-1]["Spectra"] = data_arr[:,dcp["spectra_num"]]/data_arr[:,dcp["spectra_denom"]]
            if not len(data_sets)==0:
                if len(data_sets[-1]) < 7:
                    data_sets.pop(-1)
                elif data_sets[-1]["Energy"].size != 801: # Discards uncompleted scans
                    data_sets.pop(-1)
            data_sets.append({})
            append_new = False
            data = []
        if len(line) == 0:
            append_new = True
            continue
        # Reading metadata
        elif line[0] == "#S":
            data_sets[-1]["Scan"] = filename + '-' + line[dcp["scan_name"]]
        elif line[0] == "#P4":
            data_sets[-1]["Temperature"] = int(float(line[dcp["temperature"]]))
            data_sets[-1]["Measured_temperature"] = float(line[dcp["measured temperature"]])
        elif line[0] == "#P3":
            if np.isclose(float(line[dcp["pol angle"]]), 0, atol=1e-4, rtol=0) or np.isclose(abs(float(line[dcp["pol angle"]])), np.pi, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'LH'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), np.pi/2, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'LV'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), np.pi/4, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'CP'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), -np.pi/4, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'CN'
            else:
                data_sets[-1]["Polarization"] = "npi" # This means "ni puta idea"
        elif line[0] == "#P9": # Magnetic field on or off, simple
            if np.isclose(float(line[dcp["magnet"]]), 0, atol=1e-2, rtol=0):
                data_sets[-1]["B"] = False
            else:
                data_sets[-1]["B"] = True
        elif line[0] == "#P18":
            data_sets[-1]["Incidence angle"] = float(line[dcp["incidence angle"]])
            data_sets[-1]["zSample"] = float(line[dcp["zsample"]])
        else:
            try:
                _ = float(line[0]) # If line doesn't stars with a number there is no data, it's metadata
                data.append(line)
            except ValueError:
                continue
    
    if len(data_sets[-1]) < 9: # Discards scans with non-readed metadata
        data_sets.pop(-1)
    elif data_sets[-1]["Energy"].size != 801: # Discards uncompleted scans
        data_sets.pop(-1)
    
    fl.close()
    
    # Code to discriminate samples
    # Only one strategy implemented
    
    if not sample_list:
        for d in data_sets:
            d["Sample"] = "no-name"
    
    else:
        if not disc_mode in ["region"]:
            raise ValueError("Non-valid sample discrimination mode")
            return
        if disc_mode == "region":
            if len(region_limits)+1 != len(sample_list):
                raise ValueError("Number of regions does not fit with number of samples")
                return
            region_limits = np.sort(region_limits)
            for d in data_sets:
                sample_assigned = False
                for e in exceptions:
                    try:
                        if np.isclose(e[1], d["zSample"], atol=1, rtol=0):
                            d["Sample"] = sample_list[e[0]]
                            sample_assigned = True
                    except:
                        pass
                for i,lim in enumerate(region_limits):
                    if d["zSample"] <= lim and not sample_assigned:
                        d["Sample"] = sample_list[i]
                        sample_assigned = True
                if d["zSample"] > region_limits[-1] and not sample_assigned:
                    d["Sample"] = sample_list[-1]
                    sample_assigned = True
                if not sample_assigned:
                    d["Sample"] = "no-name"
    
    return data_sets




def data_to_csv(data_file_list, csvname, overwrite=False, **kwargs):
    """
    Reads a .dat file from boreas BL29 and convert the usefull metadata and cooked data into a .csv file.

    Parameters
    ----------
    data_file_list : list(string)
        List of .dat files to cook.
    csvname : string
        Name of .csv file when cooked data is served.
    overwrite : Bool, optional
        If True overwrites the scans already written. The default is False.
    disc_mode : string, optional
        Strategy for sample distinguishment from position. The default is "region".
    sample_list : list(string), optional
        Names of the samples. If list is empty there will not be any discrimination and all samples will recive the name "no-name". The default is [].
    region_limits : list(int,float), optional
        Limits of the regions occupied by samples with temperature changes. There must be one less limit than samples. The default is [].

    """
    data_sets = []
    for flnm in data_file_list:
        data_sets += _read_boreas_file(flnm, **kwargs)
    print("Data sets generated")
    flname = csvname + ".csv"
    flexists = flname in os.listdir()
    mode = 'r+' if (flexists and not overwrite) else 'w'
    with open(flname, mode=mode, newline='') as csvfl:
        writer = csv.writer(csvfl)
        scans_readed = []
        if mode=='w':
            writer.writerow(["Scan",
                             "Temperature",
                             "Measured temperature",
                             "Polarization",
                             "B",
                             "Incidence angle",
                             "zSample",
                             "Sample",
                             "Energy"] + [" "]*800 + ["Spectra"] + [" "]*800)
        
        elif mode=='r+':
            reader = csv.reader(csvfl)
            for l,line in enumerate(reader):
                if l!=0:
                    scans_readed.append(line[0])
            
        for dset in data_sets:
            if not dset["Scan"] in scans_readed:
                writer.writerow([dset["Scan"],
                                 dset["Temperature"],
                                 dset["Measured_temperature"],
                                 dset["Polarization"],
                                 "On" if dset["B"] else "Off",
                                 dset["Incidence angle"],
                                 dset["zSample"],
                                 dset["Sample"]] + list(dset["Energy"]) + list(dset["Spectra"]))
    
    print("Data saved")
    return
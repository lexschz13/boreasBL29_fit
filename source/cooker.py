import numpy as np
import csv
import os
from .defaults import data_cook_position as dcp





def _read_boreas_file(filename):
    data_sets = []
    append_new = True
    fl = open(filename + ".dat")
    lines = [i.strip().split() for i in fl.readlines()]
    data = []
    for line in lines:
        if append_new:
            data_arr = np.array(data).astype(np.float64)
            if data_arr.ndim >= 2:
                data_sets[-1]["Energy"] = data_arr[:,dcp["energy"]]
                data_sets[-1]["Spectra"] = data_arr[:,dcp["spectra_num"]]/data_arr[:,dcp["spectra_denom"]]
            if not len(data_sets)==0:
                if len(data_sets[-1]) < 7:
                    data_sets.pop(-1)
                elif data_sets[-1]["Energy"].size != 801:
                    data_sets.pop(-1)
            data_sets.append({})
            append_new = False
            data = []
        if len(line) == 0:
            append_new = True
            continue
        elif line[0] == "#S":
            data_sets[-1]["Scan"] = filename + '-' + line[dcp["scan_name"]]
        elif line[0] == "#P4":
            data_sets[-1]["Temperature"] = int(float(line[dcp["temperature"]]))
            data_sets[-1]["Measured_temperature"] = float(line[dcp["measured temperature"]])
        elif line[0] == "#P3":
            # print(float(line[1]) / np.pi)
            if np.isclose(float(line[dcp["pol angle"]]), 0, atol=1e-4, rtol=0) or np.isclose(abs(float(line[dcp["pol angle"]])), np.pi, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'LH'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), np.pi/2, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'LV'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), np.pi/4, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'CP'
            elif np.isclose(abs(float(line[dcp["pol angle"]])), -np.pi/4, atol=1e-4, rtol=0):
                data_sets[-1]["Polarization"] = 'CN'
            else:
                data_sets[-1]["Polarization"] = "npi"
            # print(float(line[1]) / np.pi, data_sets[-1]["Polarization"])
        elif line[0] == "#P9":
            if np.isclose(float(line[dcp["magnet"]]), 0, atol=1e-2, rtol=0):
                data_sets[-1]["B"] = False
            else:
                data_sets[-1]["B"] = True
        elif line[0] == "#P18":
            data_sets[-1]["Incidence angle"] = float(line[dcp["incidence angle"]])
            data_sets[-1]["zSample"] = float(line[dcp["zsample"]])
        else:
            try:
                # data_line = float(line[0])
                data.append(line)
            except ValueError:
                continue
    
    if len(data_sets[-1]) < 8:
        data_sets.pop(-1)
    elif data_sets[-1]["Energy"].size != 801:
        data_sets.pop(-1)
    
    fl.close()
    return data_sets




def data_to_csv(data_file_list, csvname, overwrite=False):
    data_sets = []
    for flnm in data_file_list:
        data_sets += _read_boreas_file(flnm)
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
                                 dset["Sample"]] + list(dset["Energy"]) + list(dset["Spectra"]))
    
    print("Data saved")
    return

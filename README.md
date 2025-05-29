### Boreas BL29 Fitting Data Package

With data_to_csv you can take some .dat created by BL29 software of the same experiment an reorganize the data at .csv file (necessary for fitting).

To discriminate the samples with z-position (zSample) you will need to stablish the regions where those samples are moving due to the temperature changes. I recomend a first execution without assign sample names to see the positions (all scans will assign their sample name as "no-name", original, eh?). After that you will be able to identify the position ranges and stablish their limits.

With single_fit you can make a fitting indicating the name of the scan.

The function compute_gap is computing the gap for a sample, polarization, magnetic field "On" or "Off" and incidence angle.

For fittings is necessary to indicate fitting model (or mode): lorentz, gauss, voigt or general_gauss.

Do not forget to indicate the energy window for fitting.

Do not forget to edit default values to change things like initial guess or fitting bounds for a particular fitting mode.

# TOMAS-InverseModel
Objective: implement a parameter estimation routine in the TwO-Moment Aerosol Sectional (TOMAS) microphysical algorithm with particle size distribution measurements. The directory structure and contents are:

  # Code/ 
  all FORTRAN code needed to compile the model, which includes the online estimation routine. Inputs to the model and scenarios to run are defined in array.f and box.f is the main program
  
  # PreProcessing/ 
  the MATLAB script that calculates the measured "inventory" variables from the size distribution measurements used in the inverse technique. As an example, I included the raw size distribution data from Melpitz, Germany. To generate the measurement files, call Preprocessing_measurements('Melpitz').
  
  # PreProcessing/DMPS_noisy 
  the MATLAB code that calculates the uncertainty in a size distribution measured by a differential mobility particle sizer (DMPS). This code was contributed by Yuanlong Huang and I made modifications to it so that it writes values to a file representing simulated noisy measurements given a "true" size distribution, which was predicted by TOMAS (run without the estimation technique active). To generate the uncertainty / noisy simulated data, run dmps_erun.m. As an example, I included the simulated size distribution data from one scenario (#1), so the outer loop in dmps_erun.m just needs to be adjusted to determine the noisy synthetic measurements for this scenario.

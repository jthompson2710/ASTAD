# ASTAD Algorithm
 Automated Spatiotemporal Thermal Anomaly Detection (ASTAD) Algorithm 
 (previously: ASTER Automated Statistical Spatiotemporal Thermal Anomaly Detection (AASSTAD) Algorithm) 

Contributors: James O Thompson, Claudia Corradino, Tyler Leggett, Mike Ramsey

Requires: 

          - MATLAB 2021a or newer

          - Statistics and Machine Learning Toolbox
          
          - Mapping Toolbox
          
          - Also requires all the functions in the GitHub subfolder



ALGORITHM OBJECTIVE

This code performs image preprocessing and statstical anyalyis to
indenfity subtle thermal anomalies in L1T ASTER scenes.


STEP UP INFORMATION

Make sure all the functions in the GitHub subfolder are compiled and you
will need the Statistics and Machine Learning and Mapping Toolbox
installed


Place all HDF files in one directory and then fill out the information 
about the data and volcano below. This includes:

  1. Path to directory containing the HDF files
  2. Volcano Name
  3. Longitude of Volcano
  4. Latitude of Volcano


OUTCOME

After the algorithm has run a MATLAB cell data table will be saved in 
the directoy above the directory containing all the HDF files


The table format and descriptions are below in (column order), with
each row representing one HDF file:



1: FileName - The name if the HDF file

2: Dates - The date of the file (MM/DD/YYYY)

3: Times - Time in UTC (HMMSS)

4: DayNightFlag - Day and Night Flag (local time), with 1 = Day and 
                  0 = Night

5: SceneRMap - Contains the full scene map reference information for each 
               cell

6: SceneInfo - Contains the metadata for the scene from the HDF files

7-11: SceneDN_## - Contains the full scene digital number data for each 
                   band (uint16) 

12-16: FullScene_## - Contains the full scene data (gain and offset 
                      corrected) for each band (double) 

17: ImageRmap - Contains the map reference information for each cell in 
                the subset 90 x 90 pixel area around the volcano

18-22: ImageDN_## - Contains the 90 x 90 pixel subset area digital number
                    data for each band (uint16) 

23-27: ImageScene_## - Contains the 90 x 90 pixel subset area data (gain 
                      and offset corrected) for each band (double) 

28: Gabor - The results of the Gabor processing for each image (double)

29: STThreshold - The Gabor threshold used to determine anomalies 
                  (double)

30: G_ID - Mask of anomalous pixels determined from Gabor (double)

31: G-Temp - Temperature of anomalous pixels determined from Gabor 
             (double) (Kelvin)

32: InnerAnnulusRadius - The radius in pixels of the annulus calculated 
                         determining the extent of the anomalous pixels
                         away from the volcano summit(double). This is 
                         used to determine the background temperature
                         calculation.

33: G_TAB - The above background temperature of the anomalous pixels in
            Kelvin (double)

34: G_Background - The location and temperature of the background area 
                   used in the processing

35: AnomalyCategory - The anomaly category determining the anomaly status
                      for each scene, with: 
                      0 = Summit off scene, 
                      1 = Cloudy scene
                      2 = No anomaly detected
                      3 = Anomaly detected



# ASTAD Change Detection Analysis
We have included an analytical code (ASTAD Change Detection) that begins to test detection algoithms pre-, syn-, and post- eruptions by evaluating the results of the ASTAD algorithm.

METHOD OVERVIEW
1. Use the moving average to determine the underlying trend
2. Get the moving average of the signal
3. Then, to highligth increasing trend, we compute its derivative
4. Finally, since we don't look for sudden instantaneous variation but for persistent ones, i.e. changes kept for longer time, we compute its moving average again to look if the derivative is kept high for at least 5 samples
 
The threshold considered at this stage are fixed threshold and statistical. 


# 
Last updated 06/16/2023

Copyright 2023 -- James O Thompson (University of Pittsburgh)

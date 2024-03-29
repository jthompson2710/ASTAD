# Data Table Availability

This folder contains the metadata for some example data results. The actually data can be located at: http://d-scholarship.pitt.edu/id/eprint/44804

These were developed for the following manuscript in Remote Sensing of the Environment: "this will be added when published". 

This work was supported by NASA grants 80NSSC18K1001 and 80NSSC20K1336.

The data are in MATLAB .mat format as Data Tables, and can be opened in MATLAB using the load function.
The data are also available as .DAT files and can be opened in other programming/GIS software.

The data tables are in the following structure:

1: FileName - The name if the HDF file

2: Dates - The date of the file (MM/DD/YYYY)

3: Times - Time in UTC (HMMSS)

4: DayNightFlag - Day and Night Flag (local time), with 1 = Day and 0 = Night

5: R - Contains the full scene map reference information for each cell

6: Info - Contains the metadata for the scene from the HDF files

7-11: Image_## - Contains the 90 x 90 pixel subset area data (gain and offset corrected) for each band (double)

12-16: Image##_DN - Contains the 90 x 90 pixel subset area digital number data for each band (uint16)

17: Gabor - The results of the Gabor processing for each image (double)

18: STThreshold - The Gabor threshold used to determine anomalies (double)

19: G_ID - Mask of anomalous pixels determined from Gabor (double)

20: G-Temp - Temperature of anomalous pixels determined from Gabor (double) (Kelvin)

21: InnerAnnulusRadius - The radius in pixels of the annulus calculated determining the extent of the anomalous pixels away from the volcano summit(double). This is used to determine the background temperature calculation.

22: G_TAB - The above background temperature of the anomalous pixels in Kelvin (double)

23: G_Background - The location and temperature of the background area used in the processing

24: AnomalyCategory - The anomaly category determining the anomaly status for each scene, with: 0 = Summit off scene, 1 = Cloudy scene 2 = No anomaly detected 3 = Anomaly detected

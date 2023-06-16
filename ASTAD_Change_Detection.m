%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2023
% Developed by: Claudia Corradino (INGV)
% Edited by: James O Thompson (University of Pittsburgh)
%
% Automated Spatiotemporal Thermal Anomaly Detection (ASTAD) Algorithm
% Change Detection Analysis
%
% Requires: - MATLAB 2021a or newer
%           - Statistics and Machine Learning Toolbox
%           - Mapping Toolbox
%           - Also requires all the functions in the GitHub subfolder
%
% Last updated 06/16/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM OBJECTIVE
% This code is an attempt to begin to test change detection algoithms pre-,
% syn-, and post- eruptions by analysing the results of the ASTAD algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% METHOD OVERVIEW
% 1. Use the moving average to determine the underlying trend
% 2. Get the moving average of the signal
% 3. Then, to highligth increasing trend, we compute its derivative
% 4. Finally, since we don't look for sudden instantaneous variation but 
%    for persistent ones, i.e. changes kept for longer time, we compute its 
%    moving average again to look if the derivative is kept high for at 
%    least 5 samples
%
% The threshold considered at this stage are fixed threshold and 
% statistical one. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%           y_max_TAB : this is cell G_TAB in the data tables
%           time vectors : this is cell Date and Time in the data tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Determine signal         
base=y_max_TAB;
%Moving average of the signal
base2=movmean(base,[days(120),0],'omitnan','SamplePoints',time);
%Compute its derivative
base3=movmean(diff(base2),[5,0]);
%Concentrate all data into one variable
Precursor.TAB=base;
Precursor.TABsmoothed=base2;
Precursor.TABderivative=base3;
%Compute its moving average again to look if the derivative is kept high 
%for at least 5 samples
Precursor.TABIncreasedDetection_Fixed_Threshold=[base3>1.1,0];
Precursor.TABIncreasedDetection_Stat_Threshold=[base3>quantile(base3,0.85),0];

%Extract dates from data tables
Dates = Final_Gabor_Test.Date(:); 
t1 = Dates(1);
t2 = Dates(end);
t = t1:t2;
tNum= datenum(t);
DatesNum= datenum(Dates);
sz = size(Final_Gabor_Test);
for i = 1:sz(1)
    [~,ind1] = min(abs(DatesNum(i) - tNum));
    ind(i) = ind1;
end

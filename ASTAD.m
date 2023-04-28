%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2023 -- James O Thompson (University of Pittsburgh)
% Contributors: Claudia Corradino, Tyler Leggett, Mike Ramsey
%
% Automated Spatiotemporal Thermal Anomaly Detection (ASTAD) Algorithm
%
% Requires: - MATLAB 2021a or newer
%           - Statistics and Machine Learning Toolbox
%           - Mapping Toolbox
%           - Also requires all the functions in the GitHub subfolder
%
% Last updated 04/27/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM OBJECTIVE
% This code performs image preprocessing and statstical anyalyis to
% indenfity subtle thermal anomalies in L1T ASTER scenes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP UP INFORMATION
% 
% Make sure all the functions in the GitHub subfolder are compiled and you
% will need the Statistics and Machine Learning and Mapping Toolbox
% installed
%
% Place all HDF files in one directory and then fill out the information 
% about the data and volcano below. This includes:
%   1. Path to directory containing the HDF files
%   2. Volcano Name
%   3. Longitude of Volcano
%   4. Latitude of Volcano
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTCOME
%
% After the algorithm has run a MATLAB cell data table will be saved in 
% the directoy above the directory containing all the HDF files
%
% The table format and descriptions are below in (column order), with
% each row representing one HDF file:
%
% 1: FileName - The name if the HDF file
% 2: Dates - The date of the file (MM/DD/YYYY)
% 3: Times - Time in UTC (HMMSS)
% 4: DayNightFlag - Day and Night Flag (local time), with 1 = Day and 
%                   0 = Night
% 5: SceneRMap - Contains the full scene map reference information for each 
%                cell
% 6: SceneInfo - Contains the metadata for the scene from the HDF files
% 7-11: SceneDN_## - Contains the full scene digital number data for each 
%                    band (uint16) 
% 12-16: FullScene_## - Contains the full scene data (gain and offset 
%                       corrected) for each band (double) 
% 17: ImageRmap - Contains the map reference information for each cell in 
%                 the subset 90 x 90 pixel area around the volcano
% 18-22: ImageDN_## - Contains the 90 x 90 pixel subset area digital number
%                     data for each band (uint16) 
% 23-27: ImageScene_## - Contains the 90 x 90 pixel subset area data (gain 
%                       and offset corrected) for each band (double) 
% 28: Gabor - The results of the Gabor processing for each image (double)
% 29: STThreshold - The Gabor threshold used to determine anomalies 
%                   (double)
% 30: G_ID - Mask of anomalous pixels determined from Gabor (double)
% 31: G-Temp - Temperature of anomalous pixels determined from Gabor 
%              (double) (Kelvin)
% 32: InnerAnnulusRadius - The radius in pixels of the annulus calculated 
%                          determining the extent of the anomalous pixels
%                          away from the volcano summit(double). This is 
%                          used to determine the background temperature
%                          calculation.
% 33: G_TAB - The above background temperature of the anomalous pixels in
%             Kelvin (double)
% 34: G_Background - The location and temperature of the background area 
%                    used in the processing
% 35: AnomalyCategory - The anomaly category determining the anomaly status
%                       for each scene, with: 
%                       0 = Summit off scene, 
%                       1 = Cloudy scene
%                       2 = No anomaly detected
%                       3 = Anomaly detected
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc,close all
% WELCOME
uiwait(msgbox({'This code performs image preprocessing and statstical anyalyis';...
    'to identify subtle thermal anomalies in L1T ASTER scenes';'';...
    'Thank you for using our model and we hope you find it helful';'';...
    '“All models are wrong, but some are useful” George E. P. Box';''}...
    , 'Automated Spatiotemporal Thermal Anomaly Detection (ASTAD) Algorithm'));

uiwait(msgbox({'Place all HDF files in one directory and then fill out the ';...
    'information about the data and volcano below. This includes:';'';...
    '1. Path to directory containing the HDF files';...
    '2. Volcano Name';...
    '3. Longitude of Volcano';...
    '4. Latitude of Volcano';''}, 'STEP UP INFORMATION'));    
    
%% Change these variables for your dataset/volcano

runformat = 1; %This is the runnig format. Change to 0 if want to add 
%               numbers in the editor below; 1 for interactive dialog boxes

if runformat == 0
% Path to directory containing the HDF files
HDFDIR = 'D:\Users\James\OneDrive - University of Pittsburgh\MATLAB ANOMALY DETECTION CODE 1\Example Images\HDF1\';
% Volcano Name
VolcanoName = 'Etna';
% Longitude of Volcano
VolcanoLongitude = 14.9959;
% Latitude of Volcano
VolcanoLatitude = 37.7490;
end

if runformat == 1
% Directory of  HFD files
    HDFDIR = uigetdir('', 'Select directory containing the HDF files');
    HDFDIR  = join([HDFDIR '\'],'');
% DIALOG INPUT VALUES
    prompt = {'Volcano Name:','Volcano Longitude:','VolcanoLatitude:'};
    dlgtitle = 'Input Paramenters';
    definput = {'Etna','14.9959','37.7490'};
    opts.Interpreter = 'tex';
    dims = [1 65];
    input_answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
    %%assign user input values to variables
    VolcanoName = input_answer{1};
    VolcanoLongitude = str2double(input_answer{2});
    VolcanoLatitude = str2double(input_answer{3});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  It is not recommended to change anything after this line  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm Part 1: The section saves the HDF data into a table

% Saves all the images from the folder named HDF1 as a MATLAB variable 
HDFfiles = dir(join([HDFDIR '*.hdf'],'')); 
%define output directory
OUTDIR = strsplit(HDFDIR,'\');
OUTDIR = string(join(OUTDIR(1:end-2),'\'));
OUTDIR = join([OUTDIR '\'],'');

% The constants and center wavelengths to calculate Brightness Temperature 
% for each ASTER band 10-14.
k = 0;
constant1 = 119104200;
constant2 = 14387.75;
Wave10 = 8.3;
Wave11 = 8.65;
Wave12 = 9.1;
Wave13 = 10.6;
Wave14 = 11.3;

% loop through all HDF files
for g = 1:length(HDFfiles)
        
        % get each scene metadata
        iSceneInfo = hdfinfo([HDFDIR HDFfiles(g).name]);
        % determine name of TIR data
        iTIRindex = find((strcmp({iSceneInfo.Vgroup.Name}, 'TIR')),1);
        % read bands 10 to 14 
        iSceneDN_10 = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(2).SDS(1).Name);
        iSceneDN_11 = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(2).SDS(2).Name);
        iSceneDN_12 = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(2).SDS(3).Name);
        iSceneDN_13 = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(2).SDS(4).Name);
        iSceneDN_14 = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(2).SDS(5).Name);
        
        % convert band digital number to actual values using gain and
        % offset values - gives values in Celsius
        % Band 10
        iFullScene_10 = double(iSceneDN_10);
        iFullScene_10(iFullScene_10 ==0) = NaN;
        iFullScene_10 = (iFullScene_10-1)*0.006822;
        c1Wave = constant1/Wave10^5;
        iFullScene_10 = c1Wave./iFullScene_10;
        iFullScene_10 = iFullScene_10 + 1;
        iFullScene_10 = log(iFullScene_10);
        iFullScene_10 = iFullScene_10*Wave10;
        iFullScene_10 = constant2./iFullScene_10 - 273;

        % Band 11
        iFullScene_11 = double(iSceneDN_11);
        iFullScene_11(iFullScene_11 ==0) = NaN;
        iFullScene_11 = (iFullScene_11-1)*0.006780;
        c1Wave11 = constant1/Wave11^5;
        iFullScene_11 = c1Wave11./iFullScene_11;
        iFullScene_11 = iFullScene_11 + 1;
        iFullScene_11 = log(iFullScene_11);
        iFullScene_11 = iFullScene_11*Wave11;
        iFullScene_11 = constant2./iFullScene_11 - 273;

        % Band 12
        iFullScene_12 = double(iSceneDN_12);
        iFullScene_12(iFullScene_12 ==0) = NaN;
        iFullScene_12 = (iFullScene_12-1)*0.006590;
        c1Wave12 = constant1/Wave12^5;
        iFullScene_12 = c1Wave12./iFullScene_12;
        iFullScene_12 = iFullScene_12 + 1;
        iFullScene_12 = log(iFullScene_12);
        iFullScene_12 = iFullScene_12*Wave12;
        iFullScene_12 = constant2./iFullScene_12 - 273;

        % Band 13
        iFullScene_13 = double(iSceneDN_13);
        iFullScene_13(iFullScene_13 ==0) = NaN;
        iFullScene_13 = (iFullScene_13-1)*0.005693;
        c1Wave13 = constant1/Wave13^5;
        iFullScene_13 = c1Wave13./iFullScene_13;
        iFullScene_13 = iFullScene_13 + 1;
        iFullScene_13 = log(iFullScene_13);
        iFullScene_13 = iFullScene_13*Wave13;
        iFullScene_13 = constant2./iFullScene_13 - 273;

        % Band 14
        iFullScene_14 = double(iSceneDN_14);
        iFullScene_14(iFullScene_14 ==0) = NaN;
        iFullScene_14 = (iFullScene_14-1)*0.005225;
        c1Wave14 = constant1/Wave14^5;
        iFullScene_14 = c1Wave14./iFullScene_14;
        iFullScene_14 = iFullScene_14 + 1;
        iFullScene_14 = log(iFullScene_14);
        iFullScene_14 = iFullScene_14*Wave14;
        iFullScene_14 = constant2./iFullScene_14 - 273;
        
        %determine raster size of one band
        rs = size(iFullScene_10);
        %combine all bands together
        TIR = cat(3,iFullScene_10,iFullScene_11,iFullScene_12,iFullScene_13,iFullScene_14);

        % Open and the lat and long matrix associated with each scene
        GeodeticLatitude = hdfread([HDFDIR HDFfiles(g).name],  iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(1).SDS(1).Name);
        Longitude = hdfread([HDFDIR HDFfiles(g).name], iSceneInfo.Vgroup(iTIRindex).Vgroup(1).Vgroup(1).SDS(2).Name);
        
        % Save the upper left, upper right, lower left and lower right
        % lat/long coordinates into one variable inputpoints. 
        Geolat = GeodeticLatitude([1 11 111 121]);
        Geolong = Longitude([1 11 111 121]);
        inputpoints = [Geolat' Geolong'];
    
        % Save the min and max lat/long for each scene. 
        LatMax = max(inputpoints(:,1));
        LatMin = min(inputpoints(:,1));
        LonMax = max(inputpoints(:,2));
        LonMin = min(inputpoints(:,2));
    
        % Convert the lat/long to UTM
        [UTM_X_MAX,UTM_Y_MAX,~] = deg2utm(LatMax,LonMax);
        [UTM_X_MIN,UTM_Y_MIN] = deg2utm(LatMin,LonMin);
    
        % Use the min and max UTM as limits for each scene
        xWorldLimits = [UTM_X_MIN UTM_X_MAX];
        yWorldLimits = [UTM_Y_MIN UTM_Y_MAX];
        rasterSize = size(iFullScene_10);
        % Extent = 90;

        % Create the georeference information.
        try
            iSceneR = maprefcells(xWorldLimits,yWorldLimits,rasterSize,'ColumnsStartFrom','north');
             
        catch
                continue
        end
        
        % determine filename of scene 
        ifilename = erase(HDFfiles(g).name,".hdf");
        % determine data of scene 
        iDate = ifilename(12:19);
        iDate = datestr(datenum(num2str(iDate, '%d'), 'mmddyyyy'), 'mm/dd/yyyy');
        % determine UTC time of scene 
        iUTCTime = ifilename(20:25);
        iUTCTime = str2num(iUTCTime);
        % determine time zoe of scene
        zd = timezone(VolcanoLongitude);
        % determine local time of scene         
        ilocalTime = iUTCTime - zd*(10000);
        if (ilocalTime < 0); ilocalTime = ilocalTime+240000; end;
        if (ilocalTime > 240000); ilocalTime = ilocalTime-240000; end;

        % determine if scence is day or night based on local time
        iDayOrNight = [];
        if ilocalTime < 120000
                iDayOrNight = 1; %"Day";
            else
                iDayOrNight = 0; %"Night";
        end

    % move index on by one and collate all variables for each scence
    k = k+1;
    FileName{k} = ifilename;
    Dates{k} = iDate;
    Times{k} = iUTCTime;
    DayNightFlag{k} = iDayOrNight;
    SceneRmap{k} = iSceneR;
    SceneInfo{k} = iSceneInfo;
    SceneDN_10{k} = iSceneDN_10;
    SceneDN_11{k} = iSceneDN_11;
    SceneDN_12{k} = iSceneDN_12;
    SceneDN_13{k} = iSceneDN_13;
    SceneDN_14{k} = iSceneDN_14;
    FullScene_10{k} = iFullScene_10;
    FullScene_11{k} = iFullScene_11;
    FullScene_12{k} = iFullScene_12;
    FullScene_13{k} = iFullScene_13;
    FullScene_14{k} = iFullScene_14;
end

% combine all variables and convert to mat table, save out
Data1 = [FileName',Dates',Times',DayNightFlag',SceneRmap',SceneInfo',SceneDN_10',SceneDN_11',SceneDN_12',SceneDN_13',SceneDN_14',FullScene_10',FullScene_11',FullScene_12',FullScene_13',FullScene_14'];
DataTable1 = cell2table(Data1,'VariableNames',{'FileName' 'Dates' 'Times' 'DayNightFlag' 'SceneRmap' 'SceneInfo' 'SceneDN_10' 'SceneDN_11' 'SceneDN_12' 'SceneDN_13' 'SceneDN_14' 'FullScene_10' 'FullScene_11' 'FullScene_12' 'FullScene_13' 'FullScene_14'});
save(join([OUTDIR VolcanoName,'DataTableAll.mat'],''), 'DataTable1', '-v7.3')

fprintf('Finished compiling data in to a table \n');
% clear all variables except those required later
clearvars -except DataTable1 HDFfiles HDFDIR OUTDIR VolcanoName VolcanoLongitude VolcanoLatitude k constant1 constant2 Wave10 Wave11 Wave12 Wave13 Wave14    

%% Algorithm Part 2: This section crops around the target volcano

% VolcanoDiamter sets the pixel radius for the volcanoes area. 
VolcanoDiameter = 90;  
% x,y are the UTM coordinates for northing and easting, respectively.
[VolcanoX,VolcanoY,Volcanoutmzone] = deg2utm(VolcanoLatitude,VolcanoLongitude);
%Etna:499422, Popo: 539883, Fuego: 728471, Lascar: 630319, Klyuchevskoy: 602559, At: 695602
%Etna:4178154, Popo: 2103382, Fuego: 1601323, Lascar: 7415532, Klyuchevskoy: 6214671, At: 1612855

%loop through each scene -  now rows in table
for g = 1:size(DataTable1,1)
    % get filename, metadata, and map georefernece data
    ifilename = DataTable1.FileName{g};
    iSceneInfo = DataTable1.SceneInfo(g);
    iSceneRmap = DataTable1.SceneRmap(g);
    % Create a NaN buffer around the image to help perserve images where to
    % summit is located on the edge of the scene. 
    iExtentedScene_14 = DataTable1.FullScene_14{g};
        rasterSize = size(iExtentedScene_14);
        top = NaN(VolcanoDiameter,rasterSize(1,2));
        side = NaN(rasterSize(1,1)+180,VolcanoDiameter);
        iExtentedScene_14 = [top;iExtentedScene_14;top];
        iExtentedScene_14 = [side iExtentedScene_14 side];
        
    iExtentedScene_13 = DataTable1.FullScene_13{g};
        iExtentedScene_13 = [top;iExtentedScene_13;top];
        iExtentedScene_13 = [side iExtentedScene_13 side];
        
    iExtentedScene_12 = DataTable1.FullScene_12{g};
        iExtentedScene_12 = [top;iExtentedScene_12;top];
        iExtentedScene_12 = [side iExtentedScene_12 side];
        
    iExtentedScene_11 = DataTable1.FullScene_11{g};
        iExtentedScene_11 = [top;iExtentedScene_11;top];
        iExtentedScene_11 = [side iExtentedScene_11 side];
        
    iExtentedScene_10 = DataTable1.FullScene_10{g};
        iExtentedScene_10 = [top;iExtentedScene_10;top];
        iExtentedScene_10 = [side iExtentedScene_10 side];
    
    iExtentedSceneDN_14 = DataTable1.SceneDN_14{g};
        iExtentedSceneDN_14 = [top;iExtentedSceneDN_14;top];
        iExtentedSceneDN_14 = [side iExtentedSceneDN_14 side];
        
    iExtentedSceneDN_13 = DataTable1.SceneDN_13{g};
        iExtentedSceneDN_13 = [top;iExtentedSceneDN_13;top];
        iExtentedSceneDN_13 = [side iExtentedSceneDN_13 side];
        
    iExtentedSceneDN_12 = DataTable1.SceneDN_12{g};
        iExtentedSceneDN_12 = [top;iExtentedSceneDN_12;top];
        iExtentedSceneDN_12 = [side iExtentedSceneDN_12 side];
        
    iExtentedSceneDN_11 = DataTable1.SceneDN_11{g};
        iExtentedSceneDN_11 = [top;iExtentedSceneDN_11;top];
        iExtentedSceneDN_11 = [side iExtentedSceneDN_11 side];
        
    iExtentedSceneDN_10 = DataTable1.SceneDN_10{g};
        iExtentedSceneDN_10 = [top;iExtentedSceneDN_10;top];
        iExtentedSceneDN_10 = [side iExtentedSceneDN_10 side];

    % calculate the geographical corners of the extended scene and assign
    % each pixel a geolocation
    iRefMatrix = [[0 iSceneRmap.CellExtentInWorldY*-1];[iSceneRmap.CellExtentInWorldX 0];[iSceneRmap.XWorldLimits(1) iSceneRmap.YWorldLimits(2)]];
    [iX,iY] = map2pix(iRefMatrix, VolcanoX, VolcanoY);
    iX = uint16(iX);
    iY = uint16(iY);
    iX = iX + VolcanoDiameter;
    iY = iY + VolcanoDiameter;

    % Crop the ASTER L1T DN and G/O corrected scenes
    try
    iSquarecroppedImage_14 = iExtentedScene_14(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    catch
        continue
    end
    iSquarecroppedImage_13 = iExtentedScene_13(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImage_12 = iExtentedScene_12(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImage_11 = iExtentedScene_11(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImage_10 = iExtentedScene_10(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iTIR_Cropped = cat(3,iSquarecroppedImage_10,iSquarecroppedImage_11,iSquarecroppedImage_12,iSquarecroppedImage_13,iSquarecroppedImage_14);

    iSquarecroppedImageDN_14 = iExtentedSceneDN_14(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImageDN_13 = iExtentedSceneDN_13(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImageDN_12 = iExtentedSceneDN_12(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImageDN_11 = iExtentedSceneDN_11(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    iSquarecroppedImageDN_10 = iExtentedSceneDN_10(iX-VolcanoDiameter:iX+VolcanoDiameter, iY-VolcanoDiameter:iY+VolcanoDiameter);
    
    % determine new lwer and upper x and y positions 
    ixLL = iX - VolcanoDiameter;
    iyLL = iY - VolcanoDiameter;
    
    ixUL = iX + VolcanoDiameter;
    iyUL = iY + VolcanoDiameter;
 
    ixLL = double(ixLL);
    iyLL = double(iyLL);
    ixUL = double(ixUL);
    iyUL = double(iyUL);
    
    % calculate the geographical corners of the cropped scene and assign
    % each pixel a geolocation
    [iXLL,iYLL] = pix2map(iRefMatrix, ixLL, iyLL);
    [iXUL,iYUL] = pix2map(iRefMatrix, ixUL, iyUL);
    
    % determine actual geoloaction of x and y corners of cropped scene
    xWorldLimits = [iXLL iXUL];
    yWorldLimits = [iYUL iYLL];
    rasterSize = size(iSquarecroppedImage_14);
    
    % Create new georeferenced information for the cropped scene for the 
    iImageR = maprefcells(xWorldLimits,yWorldLimits,rasterSize,'ColumnsStartFrom','north');

    % collate all variables for each scence
    ImageRmap{g} = iImageR;
    ImageDN_10{g} = iSquarecroppedImageDN_10;
    ImageDN_11{g} = iSquarecroppedImageDN_11;
    ImageDN_12{g} = iSquarecroppedImageDN_12;
    ImageDN_13{g} = iSquarecroppedImageDN_13;
    ImageDN_14{g} = iSquarecroppedImageDN_14;
    Image_10{g} = iSquarecroppedImage_10;
    Image_11{g} = iSquarecroppedImage_11;
    Image_12{g} = iSquarecroppedImage_12;
    Image_13{g} = iSquarecroppedImage_13;
    Image_14{g} = iSquarecroppedImage_14;
    
end

% combine all variables and combine with initial mat table, save out
Data2 = [ImageRmap',ImageDN_10',ImageDN_11',ImageDN_12',ImageDN_13',ImageDN_14',Image_10',Image_11',Image_12',Image_13',Image_14'];
DataTable2 = cell2table(Data2,'VariableNames',{'ImageRmap' 'ImageDN_10' 'ImageDN_11' 'ImageDN_12' 'ImageDN_13' 'ImageDN_14' 'Image_10' 'Image_11' 'Image_12' 'Image_13' 'Image_14'});
DataTableAll = [DataTable1 DataTable2];
save(join([OUTDIR VolcanoName,'DataTableAll.mat'],''), 'DataTableAll', '-v7.3')

fprintf('Finished creating volcano subsets \n');
% clear all variables except those required later
clearvars -except DataTableAll VolcanoName OUTDIR

%% Algorithm Part 3: This section determines anomalies
% % %*************Main processing**********% % %
% This first runs the annulus creater to determine the background area 
% Then this runs the image preprocessing, and the Gabor filter for each cropped
% ASTER L1T scene.

% A3 is the final Gabor filter image.
% yTemps is the binary image which identifies the thermal anomalies.
% WSA1 is the temperature of the identified anomalies.
% WSA2 is the temperature above background for the identified anomalies.
% background1 is the background temperature for each scene. 

% determine table size
sz = size(DataTableAll);
TSize = size(DataTableAll);
% Sort the ASTER images in chronological order. 
DataTableAll.Dates = datenum(DataTableAll.Dates(1:sz(1)));
DataTableAll = sortrows(DataTableAll,'Dates','ascend');
DataTableAll.Dates = datetime(DataTableAll.Dates(1:sz(1)), 'ConvertFrom', 'datenum', 'Format', 'M/dd/yyyy');
counter = 1;

%% Algorithm Part 3.1: This section calculates the Gabor 

% loop through each scene -  now rows in table
for k = 1:TSize(1)
     
    %build circle around summit to make sure it is in scene.
    Temp13 = DataTableAll.Image_13{k};
    sz = size(Temp13); 
    CS = uint8(sz(1)/2); %this is used later in annulus section as well
    SC = [CS, CS, 8];
    SCth = 0:pi/50:2*pi;
    SCxunit = double(SC(3)) .* cos(SCth) + double(SC(1));
    SCyunit = double(SC(3)) .* sin(SCth) + double(SC(2));
    [C,iR] = meshgrid(1:sz(2),1:sz(1)); %this is used later in annulus section as well

    SCmask = inpolygon(iR,C,SCyunit,SCxunit);
    SCvalues = Temp13(SCmask==1);
    
    %if summit is not in scence then fill empty array and move to next scene
    if all(isnan(SCvalues(:)))
        Gabor{k} = (nan(sz(1)));
        AnomalyCategory{k} = 0; %no data over summit
        continue
    end
    AnomalyCategory{k} = 1; % data over summit
    %**********************************************************
    % Gabor calculation   
    A1 = DataTableAll.Image_13{k};
    A2 = A1;
    A2(isnan(A2)) = mean(A2(:),'omitnan');
    % prepro is a MATLAB function which runs the image preprocessing
    A2 = prepro(A2);
    if isempty(A2) == 1
        continue
    end
    % gabor_fun is a MATLAB function which runs the gabor filter              
    A3 = (A2) .* gabor_fun(A2,10);
    A3 = (A3) .* gabor_fun(A3,10);
    Gabor{k} = A3;

end

% combine all variables and combine with initial mat table
GaborDataTable = cell2table(Gabor','VariableNames',{'Gabor'});
DataTableAll = [DataTableAll GaborDataTable];

fprintf('Finished Gabor \n');
% clear all variables except those required later
clearvars -except DataTableAll VolcanoName OUTDIR TSize sz AnomalyCategory counter

%% Algorithm Part 3.2: The section collates Gabor values and determines threshold

% assign standard deviation of gabor detection
ns=3;
% loop through each scene -  now rows in table
for k = 1:TSize(1)
    if isempty(DataTableAll.Gabor{k})== 0
        allGabor(:,:,k) = DataTableAll.Gabor{k};
    else
        allGabor(:,:,k) = nan(size(DataTableAll.Image_13{k}));
    end
end

% determine gabor threshold value
STThreshold = mean(allGabor(:),'omitnan')+ns.*std(allGabor(:),'omitnan');
STThresholdAll = num2cell(ones(1,TSize(1)) .* STThreshold);

% combine all variables with initial mat table
STThresholdAllTable = cell2table(STThresholdAll','VariableNames',{'STThreshold'});
DataTableAll = [DataTableAll STThresholdAllTable];
DataTableAll.STThreshold = STThresholdAll';

fprintf('Finished Thresholding \n');
% clear all variables except those required later
clearvars -except DataTableAll VolcanoName OUTDIR TSize sz AnomalyCategory counter

%% Algorithm Part 3.3: The section determines anomalies based on 
%  Gabor values and the threshold

% loop through each scene -  now rows in table
for k = 1:TSize(1)
    % apply gabor threshold to gabor filter
    iGabor = DataTableAll.Gabor{k};
    if isnan(iGabor)
        ID_Points{k} = (nan(sz(1)));
        ID_Temp{k} = (nan(sz(1)));
        continue
    end
    % find pixels with values greater than threshold
    idx = find(iGabor(:) > DataTableAll.STThreshold{k});
    % assign these pixels a mask value
    yTemps = iGabor;
    yTemps(idx) = 100000000;
    yTemps(yTemps < 100000000) = 0;
    yTemps(yTemps == 100000000) = 1;
    yTemps(isnan(yTemps)) = 0;
    % replace all pixels with no anomaly with nan
    WSA1 = DataTableAll.Image_13{k};
    WSA1(yTemps == 0) = NaN;
	
    % collate pixels and temperature values of anomalies
	ID_Points{k} = yTemps;
    ID_Temp{k} = WSA1;
end

% combine all variables and combine with initial mat table
IDDataTable = cell2table([ID_Points', ID_Temp'],'VariableNames',{'G_ID' 'G_Temp'});
DataTableAll = [DataTableAll IDDataTable];

fprintf('Finished applying Gabor Thresholding to determine anomalies \n');
% clear all variables except those required later
clearvars -except DataTableAll VolcanoName OUTDIR TSize sz AnomalyCategory counter

%% Algorithm Part 3.4: The section creates annulus, calculates background,
%  and calculates the temperature of detected anomalies

% loop through each scene -  now rows in table
for k = 1:TSize(1)
	%load computed gabor above
	iGabor = DataTableAll.Gabor{k};  
    % skip if no anomalies previously detected
    if isnan(iGabor)
        ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
		ICradius{k} = nan;
        AnomalyCategory{k} = 1; % cloudy 
        counter = counter + 1;
        continue
    end
    % create intial annulus size
    sz = size(DataTableAll.Image_13{k});
	sd = zeros(1, floor((sz(1)/2-5)));
    CS = uint8(sz(1)/2); 
    [C,iR] = meshgrid(1:sz(2),1:sz(1)); 

    % skip if no anomalies previously detected
	if all(DataTableAll.G_ID{k},'all')
        ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
		ICradius{k} = nan;
        AnomalyCategory{k} = 1; % cloudy 
        counter = counter + 1;
        continue		
    end

    % loop through annulus radii to edge of cropped area
    for i = 1:((sz(1)/2)-5)
        % determine inner and outer annulus radius
        InnerRadius = i;
        OuterRadius = i+5;

        % create outer annulus
        OC = [CS, CS, OuterRadius];
        OCth = 0:pi/50:2*pi;
        OCxunit = double(OC(3)) .* cos(OCth) + double(OC(1));
        OCyunit = double(OC(3)) .* sin(OCth) + double(OC(2));
        
        % create inner annulus
        IC = [CS, CS, InnerRadius];
        ICth = 0:pi/50:2*pi;
        ICxunit = double(IC(3)) .* cos(ICth) + double(IC(1));
        ICyunit = double(IC(3)) .* sin(ICth) + double(IC(2));

        % create annulus mask
        mask = inpolygon(iR,C,OCyunit,OCxunit) & ~inpolygon(iR,C,ICyunit,ICxunit);
        croppedImage = iGabor(:,:,1).*double(mask);
        croppedImage(croppedImage == 0)= NaN; 

        % calculate annulus MAD
        av = mad(croppedImage(:),1);
        sd(i) = av;
    end
    
    sdmad=sd;
    sdmad(isnan(sdmad))=0;
    TFmad = islocalmin(sdmad); %determine local min
    sdmad87 = prctile(sdmad,87);
    xx = 1:length(sd);
    TFmadxx = xx(TFmad);
    % if no TFmadxx and so no hotspots then fill empty array and move to next scene
    if isempty(TFmadxx)
        ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
        ICradius{k} = nan;
        AnomalyCategory{k} = 2; %no anomlay detected
		counter = counter + 1;
        continue
    end
        
    % use local minimium and 87% threshold to determine annulus extent
    j=1;
    while any(sdmad(TFmadxx(j):end) > sdmad87)
        j=j+1;
        if j > length(TFmadxx)
            j=0;
            break
        end
    end
    
    if j==0
        ICradius{k} = 2; % save inner radius - this is the minimum value
    else
        ICradius{k} = TFmadxx(j); % save inner radius
    end
    

    annulusinnerradius = ICradius{k};
    % if not annulus created, skip
    if isempty(annulusinnerradius)== 1 
        ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
        AnomalyCategory{k} = 2; %no anomlay detected
		counter = counter + 1;
        continue
    end

    % use the background annulus mask conputed earlier
    % create background outer annulus (= inner annulus +5)
    backgroundOC = [CS, CS, (annulusinnerradius+5)];
    backgroundOCth = 0:pi/50:2*pi;
    backgroundOCxunit = double(backgroundOC(3)) .* cos(backgroundOCth) + double(backgroundOC(1));
    backgroundOCyunit = double(backgroundOC(3)) .* sin(backgroundOCth) + double(backgroundOC(2));
    % create background inner annulus
    backgroundIC = [CS, CS, annulusinnerradius];
    backgroundICth = 0:pi/50:2*pi;
    backgroundICxunit = double(backgroundIC(3)) .* cos(backgroundICth) + double(backgroundIC(1));
    backgroundICyunit = double(backgroundIC(3)) .* sin(backgroundICth) + double(backgroundIC(2));
    % create background annulus mask
    backgroundmask = inpolygon(iR,C,backgroundOCyunit,backgroundOCxunit) & ~inpolygon(iR,C,backgroundICyunit,backgroundICxunit);
    % extract only background values in mask
    backgroundvalues = DataTableAll.Image_13{k} .* backgroundmask;
    % skip if no background needed 
    check  = isempty(find(isnan(backgroundvalues) == 0)); 
    if check == 1
		ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
        AnomalyCategory{k} = 2; %no anomlay detected
		counter = counter + 1;
        continue
    else
    % calculate mean background value    
    backgroundvalues(backgroundvalues == 0) = NaN;
    backgroundvalues = rmoutliers(backgroundvalues,'mean');
    % calculate temperature above background of any anomalous pixels
    CalTAB = DataTableAll.Image_13{k} - mean(backgroundvalues(:),'omitnan');
    croppedTAB = CalTAB .* inpolygon(iR,C,backgroundOCyunit,backgroundOCxunit);
    DataTableAll.G_Temp{k} = DataTableAll.G_Temp{k} .* inpolygon(iR,C,backgroundOCyunit,backgroundOCxunit);
    DataTableAll.G_Temp{k}(DataTableAll.G_Temp{k} == 0) = NaN;
    DataTableAll.G_ID{k} = DataTableAll.G_ID{k} .* inpolygon(iR,C,backgroundOCyunit,backgroundOCxunit);
    DataTableAll.G_ID{k}(DataTableAll.G_ID{k} == 0) = NaN;
    % remove any negative values
    croppedTAB(croppedTAB <= 0) = NaN;
    WSA2 = croppedTAB;
    WSA2(isnan(DataTableAll.G_Temp{k})) = NaN;
    
    % assign values if no anomalies
    if all(isnan(WSA2(:)))
        ID_TAB{k} = WSA2;
        Background{k} = backgroundvalues;
        AnomalyCategory{k} = 2; %no anomlay detected
        counter = counter + 1;
        continue
    end
    
    % assign anomaly values
    ID_TAB{k} = WSA2;
    Background{k} = backgroundvalues;
    AnomalyCategory{k} = 3; % anomlay detected

    fprintf('Finished background iteration #%d\n', counter);
    counter = counter + 1;
    end

end

%% Algorithm Part 4: The section save out all variables

Outputdata = [ICradius', ID_TAB', Background', AnomalyCategory'];
Outputdatatable = cell2table(Outputdata,'VariableNames',{'InnerAnnulusRadius' 'G_TAB' 'G_Background' 'AnomalyCategory'});
DataTableAll = [DataTableAll Outputdatatable];
save(join([OUTDIR VolcanoName,'DataTableAll.mat'],''), 'DataTableAll', '-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


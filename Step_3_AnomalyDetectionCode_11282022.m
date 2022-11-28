%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 -- James O Thompson (University of Pittsburgh)
% Contributors: Claudia Corradino (INGV), Tyler Leggett (Pitt), 
% Mike Ramsey (Pitt)
%
% ASTER Automated Statistical Spatiotemporal Thermal Anomaly 
% Detection (AASSTAD) Algorithm
%
% This code performs image preprocessing and statstical anyalyis to
% indenfity subtle thermal anomalies in the cropped L1T ASTER scene.
%
%Requires MATLAB 2021a or newer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc,close all
%% Set inizialization parameters
volcano='Klyuchevskoy';
ns=2; %number of std for the temporal threshold
path=cd; %if custom, insert it here: path='D:\Users\Name\MATLAB ANOMALY DETECTION CODE'
check_gabor=0;% 1 if Gabor already computed, 0 otherwise
%% Loading data
load([path,'\Full Dataset\',volcano,'\Final_',volcano,'.mat'])
sz = size(Final);
% Sort the ASTER images in chronological order. 
Final.Date = datenum(Final.Date(1:sz(1)));
Final = sortrows(Final,'Date','ascend');
Final.Date = datetime(Final.Date(1:sz(1)), 'ConvertFrom', 'datenum', 'Format', 'M/dd/yyyy');
%% Day-Night flag
% This sets the time threshold to create the day and night flag. Night is
% giving the value of 0, and day is given 1. 
f1 = str2func('@(x,y) x>y');
f2 = str2func('@(x,y) x<y');

switch volcano
    case 'Etna'
    f=f2;
    time= 110000; 
    case 'Fuego'
    f=f1;
    time=160000; 
    case 'Klyuchevskoy'
    f=f2;
    time=100000; 
    case 'Popocatpetl'
    f=f1;
    time=100000; 
    case 'Lascar'
    f=f1;
    time=100000;
end

for i = 1:sz(1)
    if f(Final.Time{i},time)
    Final.DayNightFlag{i} = 1; %day
    else 
      Final.DayNightFlag{i} = 0; %night
    end
end
%% *************Main processing**********

% Then this runs the image preprocessing, and the Gabor filter for each cropped
% ASTER L1T scene.

% A3 is the final Gabor filter image.
% yTemps is the binary image which identifies the thermal anomalies.
% WSA1 is the temperature of the identified anomalies.
% WSA2 is the temperature above background for the identified anomalies.
% background1 is the background temperature for each scene. 

counter = 1;

tic
% figure()
if check_gabor==0
for k = 1:sz(1)
    %build circle around summit to make sure it is in scene.
    Temp13 = Final.Image_13{k};
    szT = size(Temp13); 
    CS = uint8(szT(1)/2); %this is used later in annulus section as well
    SC = [CS, CS, 8];
    SCth = 0:pi/50:2*pi;
    SCxunit = double(SC(3)) .* cos(SCth) + double(SC(1));
    SCyunit = double(SC(3)) .* sin(SCth) + double(SC(2));
    [C,R] = meshgrid(1:szT(2),1:szT(1)); %this is used later in annulus section as well

    SCmask = inpolygon(R,C,SCyunit,SCxunit);
    SCvalues = Temp13(SCmask==1);
    
    %if summit is not in scence then fill empty array and move to next scene
    if all(isnan(SCvalues(:)))
        Gabor{k} = (nan(szT(1)));
        AnomalyCategory{k} = 0; %no data over summit
        continue
    end
    
    %**********************************************************
    %the gabor stuff has now been moved here    
    A1 = Final.Image_13{k};
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
    Gabor{k} = (A3);
end


GaborDataTable = cell2table(Gabor','VariableNames',{'Gabor'});
Final = [Final GaborDataTable];

fprintf('Just finished Gabor');
end
%%
%Collate Gabor stuff and determine threshold

for k = 1:sz(1)
    if isempty(Final.Gabor{k})== 0
        allGabor(:,:,k) = Final.Gabor{k};
    else
        allGabor(:,:,k) = nan(size(Final.Image_13{k}));
    end
end

STThreshold = mean(allGabor(:),'omitnan')+ns.*std(allGabor(:),'omitnan');
Final.STThreshold(:)=STThreshold;

fprintf('Just finished Thresholding');
%%

for k = 1:sz(1)
% apply gabor threshold to gabor filter
    Gabor1 = Final.Gabor{k};
    if isnan(Gabor1)
        ID_Points{k} = (nan(sz(1)));
        ID_Temp{k} = (nan(sz(1)));
        continue
    end
    yTemps = Gabor1> Final.STThreshold(k);
    WSA1 = Final.Image_13{k};
    WSA1(yTemps == 0) = NaN;	
	ID_Points{k} = (yTemps);
    ID_Temp{k} = WSA1;

fprintf('Just finished applying Gabor Thresholding to determine anomalies');

% ---------Create Annulus and Background---------------------------------
	%load computed gabor above
    if isnan(Gabor1)
        ID_TAB{k} = (nan(sz(1)));
        Background{k} = (nan(sz(1)));
		ICradius{k} = nan;
        AnomalyCategory{k} = 1; % cloudy maybe
        counter = counter + 1;
        continue
    end
    szG = size(Gabor1); 
	sd = zeros(1, floor((szG(1)/2-5)));
    CS = uint8(szG(1)/2); %this is used later in annulus section as well
    [C,R] = meshgrid(1:szG(2),1:szG(1)); %this is used later in annulus section as well

	if all(ID_Points{k},'all')
        ID_TAB{k} = (nan(szG(1)));
        Background{k} = (nan(szG(1)));
		ICradius{k} = nan;
        AnomalyCategory{k} = 1; % cloudy maybe
        counter = counter + 1;
        continue		
	end

    for i = 1:((szG(1)/2)-5)
        InnerRadius = i;
        OuterRadius = i+5;

        %create outer annulus
        OC = [CS, CS, OuterRadius];
        OCth = 0:pi/50:2*pi;
        OCxunit = double(OC(3)) .* cos(OCth) + double(OC(1));
        OCyunit = double(OC(3)) .* sin(OCth) + double(OC(2));
        
        %create inner annulus
        IC = [CS, CS, InnerRadius];
        ICth = 0:pi/50:2*pi;
        ICxunit = double(IC(3)) .* cos(ICth) + double(IC(1));
        ICyunit = double(IC(3)) .* sin(ICth) + double(IC(2));

        %create annulus mask
        mask = inpolygon(R,C,OCyunit,OCxunit) & ~inpolygon(R,C,ICyunit,ICxunit);
        croppedImage = Gabor1.*double(mask);
        croppedImage(croppedImage == 0)= NaN; 

        %calculate annulus MAD
        av = mad(croppedImage(:),1);
        sd(i) = av;
    end
    
    sdmad=sd;
    sdmad(isnan(sdmad))=0;
    TFmad = islocalmin(sdmad); %determine local min
    %sdmad95 = prctile(sdmad,95);
    sdmad87 = prctile(sdmad,87);
    %sdmad68 = prctile(sdmad,68);
    xx = 1:length(sd);
    TFmadxx = xx(TFmad);
    %if no TFmadxx and so hotspots then fill empty array and move to next scene
    if isempty(TFmadxx)
        ID_TAB{k} = (nan(szG(1)));
        Background{k} = (nan(szG(1)));
        ICradius{k} = nan;
        AnomalyCategory{k} = 2; %no anomlay detected
		counter = counter + 1;
        continue
    end
        
    %use local minimium and 87% threshold to determine annulus extent
    j=1;
    while any(sdmad(TFmadxx(j):end) > sdmad87)
        j=j+1;
        if j > length(TFmadxx)
            j=0;
            break
        end
    end
    
    if j==0
        ICradius{k} = 2; %save inner radius - this is the minimum value
    else
        ICradius{k} = TFmadxx(j); %save inner radius
    end
    

    annulusinnerradius = ICradius{k};
    %use the background annulus mask conputed earlier
    %create background outer annulus (= inner annulus +5)
    backgroundOC = [CS, CS, (annulusinnerradius+5)];
    backgroundOCth = 0:pi/50:2*pi;
    backgroundOCxunit = double(backgroundOC(3)) .* cos(backgroundOCth) + double(backgroundOC(1));
    backgroundOCyunit = double(backgroundOC(3)) .* sin(backgroundOCth) + double(backgroundOC(2));
    %create background inner annulus
    backgroundIC = [CS, CS, annulusinnerradius];
    backgroundICth = 0:pi/50:2*pi;
    backgroundICxunit = double(backgroundIC(3)) .* cos(backgroundICth) + double(backgroundIC(1));
    backgroundICyunit = double(backgroundIC(3)) .* sin(backgroundICth) + double(backgroundIC(2));
    %create background annulus mask
    backgroundmask = inpolygon(R,C,backgroundOCyunit,backgroundOCxunit) & ~inpolygon(R,C,backgroundICyunit,backgroundICxunit);

    %the old bit....
    backgroundvalues = Final.Image_13{k} .* backgroundmask;
    check  = isempty(find(isnan(backgroundvalues) == 0)); 
    if check == 1
		ID_TAB{k} = (nan(szG(1)));
        Background{k} = (nan(szG(1)));
        AnomalyCategory{k} = 2; %no anomlay detected
		counter = counter + 1;
        continue
    else
    backgroundvalues(backgroundvalues == 0) = NaN;
    backgroundvalues = rmoutliers(backgroundvalues,'mean');
    CalTAB = Final.Image_13{k} - mean(backgroundvalues(:),'omitnan');
    croppedTAB = CalTAB .* inpolygon(R,C,backgroundOCyunit,backgroundOCxunit);
    croppedTAB(croppedTAB <= 0) = NaN;
    WSA2 = croppedTAB;
    WSA2(isnan(ID_Temp{k})) = NaN;
    
    if all(isnan(WSA2(:)))
        ID_TAB{k} = WSA2;
        Background{k} = backgroundvalues;
        AnomalyCategory{k} = 2; %no anomlay detected
        counter = counter + 1;
        continue
    end
    
    ID_TAB{k} = WSA2;
    Background{k} = backgroundvalues;
    AnomalyCategory{k} = 3; % anomlay detected

    fprintf('Just finished background iteration #%d\n', counter);
    counter = counter + 1;
    end

    
    
    %**********************************************************

end
toc
%%
% save out all variables
Outputdata = [ID_Points', ID_Temp',ICradius', ID_TAB', Background', AnomalyCategory'];
Outputdatatable = cell2table(Outputdata,'VariableNames',{'G_ID' 'G_Temp','InnerAnnulusRadius' 'G_TAB' 'G_Background' 'AnomalyCategory'});
Final_Gabor_Test = [Final Outputdatatable];
save('Final_Gabor_Etna.mat', 'Final_Gabor', '-v7.3')

%% ***************                              ****************** %%

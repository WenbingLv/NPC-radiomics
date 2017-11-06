function [SUVmax,SUVpeak,SUVmean,MATV,TLG] = getSUVfeatures(ROIonlyPET,pixelW,sliceS)
% get the SUVmax, SUVpeak, SUVmean,metabolically active tumor volume (MATV) 
% and total lesion glycolysis (TLG) of the tumor volume.
%
%ljlubme@gmail.com
%Southern Medical University
%
% Initialization
ROIonlyPET = padarray(ROIonlyPET,[1 1 1],NaN);

% SUVmax
SUVmax = max(ROIonlyPET(:));

% SUVpeak (using 26 neighbors around all voxeles in ROI)
[x,y,z] = size(ROIonlyPET);
temp = 1;
for i = 2:x-1
for j = 2:y-1
for k = 2:z-1
cube26 = ROIonlyPET(i-1:i+1,j-1:j+1,k-1:k+1);
meancube26 = mean(cube26(~isnan(cube26)));
listall(temp) = meancube26;
temp = temp+1;
end
end
end
SUVpeak = max(listall);

% SUVmean
SUVmean=mean(ROIonlyPET(~isnan(ROIonlyPET)));

% MATV
mask = ~isnan(ROIonlyPET); 
numberVoxel = sum(mask(:));
MATV = numberVoxel * pixelW * pixelW * sliceS;
MATV = MATV./1e3;%cm^3

%TLG
TLG = SUVmean.*MATV;
end

% This is a matlab demo showing how to construct radiomic matrix 
%and extract radiomic features using different strategies/parameter 
%settings for PET/CT scaned NPC images. 
%
% W.Lv, Q.Yuan et al. Impact of parameter settings for radiomic features on their robustness 
% and diagnostic task performance: application to nasopharyngeal 18F-FDG PET/CT.
% European Radiology,under review.
%
%ljlubme@gmail.com
%Southern Medical University
%

clear;clc;

%% parameter
bin = 0.1;% discretization bin width
pixelW = 0.98; sliceS = 3;% - pixelW: Pixel width, or in-plane resolution, in mm.% - sliceS: Slice spacing, in mm.
neighbor = [6,18,26];% beighborhoods used in the construction of GLSZM
window = [3,5,7,9,11];% window size used in the construction of NGTDM

%% load PET (SUV images)
load PET.mat;

%% load segmentation (1: lesion, 0: non-lesion)
load T.mat;
        
id = find(T==0);
PET(id) = NaN; 

%% GLRLM_M1     
[ROIonlyM,levelsM] = Volume_discretization(PET,T,pixelW,sliceS,'pixelW',bin);%fixed bin width
GLRLM = getGLRLM_M1(ROIonlyM,levelsM);
featuresRL = getGLRLMfeatures(GLRLM);
featuresRL_M1 = featuresRL;

%% GLRLM_M13 
GLRLM = getGLRLM_M13(ROIonlyM,levelsM);
for i = 1:13 % 13 directions
featuresRL = getGLRLMfeatures(GLRLM{i});
featuresRL_M13{i} = featuresRL;
end

%% SUV features  
SUVbox = getSUVbox(PET,T);
[SUVmax,SUVpeak,SUVmean,MATV,TLG] = getSUVfeatures(SUVbox,pixelW,sliceS);

%% GLCM_1S 
[ROIonlyM,levelsM] = Volume_discretization(PET,T,pixelW,sliceS,'pixelW',bin);%fixed bin width
levelsM = max(levelsM);
idnan = find(isnan(ROIonlyM));
ROIonlyM(idnan) = levelsM + 1;
coocMatS = getGLCM_Symmetric(ROIonlyM,'distance',[1,2,3,4,5,6,7,8,9,10],'direction', [1 0 0; 1 1 0; 0 1 0; -1 1 0; 0 0 1; 1 0 1; -1 0 1; 0,1,1; 0 -1 1; 1 1 1; -1 1 1; 1 1 -1; 1 1 -1],'numgray',levelsM + 1);
for i = 1:10 %distance 
new_coMat=coocMatS(1:end-1,1:end-1,:,i);
sum_coMat=sum(new_coMat,3);
featuresC = getGLCMfeatures(sum_coMat);   
featuresC_1S{i} = featuresC; 
end

%% GLCM_1A     
coocMatA =getGLCM_Asymmetric(ROIonlyM,'distance',[1,2,3,4,5,6,7,8,9,10],'direction', [1 0 0; 1 1 0; 0 1 0; -1 1 0; 0 0 1; 1 0 1; -1 0 1; 0,1,1; 0 -1 1; 1 1 1; -1 1 1; 1 1 -1; 1 1 -1],'numgray',levelsM + 1);
for i = 1:10 %distance   
new_coMat=coocMatA(1:end-1,1:end-1,:,i);
sum_coMat=sum(new_coMat,3);% only generate a signal matrix (1A)
featuresC = getGLCMfeatures(sum_coMat);
featuresC_1A{i} = featuresC;
end

%% GLCM_13S     
for i = 1:13 %direction
    for j = 1:10 %distance          
        new_coMat = coocMatS(1:end-1,1:end-1,i,j);
        featuresC = getGLCMfeatures(new_coMat);  
        featuresC_13S{i,j} = featuresC;% need to be averaged along 13 directions for each features
    end 
end

%% GLCM_13A     
for i = 1:13 %direction
    for j = 1:10 %distance   
        new_coMat = coocMatA(1:end-1,1:end-1,i,j);
        featuresC = getGLCMfeatures(new_coMat);  
        featuresC_13A{i,j} = featuresC;% need to be averaged along 13 directions for each features
    end 
end

%% GLRLM_M1     
[ROIonlyM,levelsM] = Volume_discretization(PET,T,pixelW,sliceS,'pixelW',bin);%fixed bin width
GLRLM = getGLRLM_M1(ROIonlyM,levelsM);
featuresRL = getGLRLMfeatures(GLRLM);
featuresRL_M1 = featuresRL;

%% GLRLM_M13 
GLRLM = getGLRLM_M13(ROIonlyM,levelsM);
for i = 1:13 % 13 directions
featuresRL = getGLRLMfeatures(GLRLM{i});
featuresRL_M13{i} = featuresRL;
end

%% GLSZM_N     
for i = 1:length(neighbor)
    GLSZM = getGLSZM_neighbor(ROIonlyM,levelsM,neighbor(i));
    featuresSZ = getGLSZMfeatures(GLSZM);
    featuresSZ_N{i} = featuresSZ;
end

%% NGTDM_W
for i = 1:length(window)
    [NGTDM,countValid] = getNGTDM_window(ROIonlyM,levelsM,window(i));
    featuresTD = getNGTDMfeatures(NGTDM,countValid);
    featuresTD_W{i}= featuresTD;
end


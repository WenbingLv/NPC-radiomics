function [ROIonly,levels] = Volume_discretization(volume,mask,pixelW,sliceS,scale,bin)
% 1. Computation of the smallest box containing region of interest (ROI), 
%    if necessary (ROIbox).
% 2. Isotropic resampling.
% 3. Tumor discretization using equal bin width.
%
%ljlubme@gmail.com
%Southern Medical University
%
% COMPUTATION OF THE SMALLEST BOX CONTAINING THE ROI
[boxBound] = computeBoundingBox(mask);
maskBox = mask(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));
ROIbox = volume(boxBound(1,1):boxBound(1,2),boxBound(2,1):boxBound(2,2),boxBound(3,1):boxBound(3,2));

% ISOTROPIC RESAMPLING
flagPW = 0;
if strcmp(scale,'pixelW')
    flagPW = 1;
end
if flagPW
    a = 1;
    b = 1;
    c = sliceS/pixelW;
else
    a = pixelW/scale;
    b = pixelW/scale;
    c = sliceS/scale;
end
if numel(size(ROIbox))==3
    if a + b + c ~= 3 % If false, no resampling is needed
        maskBox = imresize3D(maskBox,[],[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b),round(double(size(maskBox,3))*c)],'nearest','fill');
        ROIbox = imresize3D(ROIbox,[],[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b),round(double(size(ROIbox,3))*c)],'cubic','fill');
    end
elseif numel(size(ROIbox))==2
    if a + b ~= 2 % If false, no resampling is needed
        maskBox = imresize(maskBox,[round(double(size(maskBox,1))*a),round(double(size(maskBox,2))*b)],'nearest');
        ROIbox = imresize(ROIbox,[round(double(size(ROIbox,1))*a),round(double(size(ROIbox,2))*b)],'cubic','Antialiasing',true);
    end
end

ROIonly = ROIbox; ROIonly(~maskBox) = NaN; ROIonly(maskBox<0) = NaN;
ROIonly(maskBox==0) = NaN;

% discretization method
[ROIonly,levels] = binWidthDiscretization(ROIonly,bin);

end

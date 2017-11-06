function [ROIonlyAbsolute,levels] = binWidthDiscretization(ROIonly,bin)
% tumor discretization using equal bin width
%
%ljlubme@gmail.com
%Southern Medical University
%
ROIonly = double(ROIonly);
minVal = min(ROIonly(:));
ROIonlyAbsolute = round((ROIonly./bin) - (minVal./bin))+1;
levels = 1:max(ROIonlyAbsolute(:));
end
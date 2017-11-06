function GLSZM = getGLSZM_neighbor(ROIOnly,levels,neighbor)

% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% - Revision: Oct 2017

% PRELIMINARY
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
levelTemp = max(levels) + 1;
ROIOnly(isnan(ROIOnly)) = levelTemp;
levels = [levels,levelTemp];


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVect = round(levels*adjust)/adjust;
ROIOnly = round(ROIOnly*adjust)/adjust;
NL = length(levels) - 1;


% INITIALIZATION
nInit = numel(ROIOnly);
GLSZM = zeros(NL,nInit);


% COMPUTATION OF GLSZM
temp = ROIOnly;
for i = 1:NL
    temp(ROIOnly~=uniqueVect(i)) = 0;
    temp(ROIOnly==uniqueVect(i)) = 1;
    connObjects = bwconncomp(temp,neighbor);
    nZone = length(connObjects.PixelIdxList);
    for j = 1:nZone
        col = length(connObjects.PixelIdxList{j});
        GLSZM(i,col) = GLSZM(i,col) + 1;
    end
end


% REMOVE UNECESSARY COLUMNS
stop = find(sum(GLSZM),1,'last');
GLSZM(:,(stop+1):end) = [];

end
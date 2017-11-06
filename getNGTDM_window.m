function [NGTDM,countValid] = getNGTDM_window(ROIOnly,levels,window)

% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% - Revision: Oct 2017

% PRELIMINARY
if numel(size(ROIOnly)) == 2 % generalization to 2D inputs
    twoD = 1;
else
    twoD = 0;
end
nLevel = length(levels);
if nLevel > 100
    adjust = 10000;
else
    adjust = 1000;
end
if twoD
    ROIOnly = padarray(ROIOnly,[1,1],NaN,'both');
else
    ROIOnly = padarray(ROIOnly,[1,1,1],NaN,'both');
end


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVol = round(levels*adjust)/adjust;
ROIOnly = round(ROIOnly*adjust)/adjust;
NL = length(levels);
temp = ROIOnly;
for i = 1:NL
    ROIOnly(temp==uniqueVol(i)) = i;
end


% INTIALIZATION
NGTDM = zeros(NL,1);
countValid = zeros(NL,1);
win = (window-1)./2;
win = int8(win);

[Sx,Sy,Sz] = size(ROIOnly);
ROItemp = zeros(Sx+2*win,Sy+2*win,Sz+2*win);
ROItemp(win+1:Sx+win,win+1:Sy+win,win+1:Sz+win) = ROIOnly;
ROIOnly = ROItemp;
id = find(ROIOnly == 0);
ROIOnly(id) = NaN;


% COMPUTATION OF NGTDM
if twoD
    [a,b] = size(ROIOnly);
    [i,j] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
    posValid = [i,j];
    nValid_temp = size(posValid,1);
    for n = 193:nValid_temp
        weights = zeros(window, window);
        if posValid(n,1) > win && posValid(n,1) < a-win+1 && posValid(n,2) > win && posValid(n,2) < b-win+1
        neighbours = zeros(window^2,1);
        neighbours(1:window^2) = ROIOnly((posValid(n,1)-win):(posValid(n,1)+win),(posValid(n,2)-win):(posValid(n,2)+win));
        centerx = posValid(n,1); centery = posValid(n,2);
        xx = 1;
       for x = (posValid(n,1)-win):(posValid(n,1)+win)
           yy = 1;
           for y = (posValid(n,2)-win):(posValid(n,2)+win) 
               x  = double(x); y = double(y); 
               if x == centerx && y == centery
               weights(xx,yy) = 1;
               else
               weights(xx,yy) = 1./sqrt((x-centerx)^2 + (y-centery)^2);
               end
               yy = yy+1;
           end
           xx = xx+1;
       end 
       weights = weights(:);
        neighbours = neighbours.*weights;
        value = neighbours((window^2+1)./2);
        neighbours((window^2+1)./2) = NaN;
        neighbours = neighbours/sum(weights(~isnan(neighbours)));
        neighbours((window^2+1)./2) = []; % Remove the center voxel
        if ~isempty(neighbours(~isnan(neighbours))) % Thus only excluding voxels with NaNs only as neighbors.
            NGTDM(value) = NGTDM(value) + abs(value-sum(neighbours(~isnan(neighbours))));
            countValid(value) = countValid(value) + 1;
        end
        else
        end
    end
else
    [i,j,k] = ind2sub(size(ROIOnly),find(~isnan(ROIOnly)));
    posValid = [i,j,k];
    nValid_temp = size(posValid,1);
    %% comput the weights for different window size
    weights = zeros(window, window, window);
    center = win + 1;
    center = double(center);
    xx = 1;
    for x = center-win:center+win
        yy = 1;
        for y = center-win:center+win
            zz = 1;
            for z = center-win:center+win
                x = double(x);y = double(y);z = double(z);
               if x == center && y == center && z == center
               weights(xx,yy,zz) = 1;
               else
               weights(xx,yy,zz) = 1./sqrt((x-center)^2 + (y-center)^2 + (z-center)^2);
               end
               zz = zz+1;
            end 
            yy = yy+1;
        end
        xx = xx+1;
    end
    weights = weights(:);   
    for n = 1:nValid_temp
        neighbours = zeros(window^3,1);
        neighbours(1:window^3) = ROIOnly((posValid(n,1)-win):(posValid(n,1)+win),(posValid(n,2)-win):(posValid(n,2)+win),(posValid(n,3)-win):(posValid(n,3)+win));
        neighbours = neighbours.*weights;
        value = neighbours((window^3-1)./2+1);
        neighbours((window^3-1)./2+1) = NaN;
        neighbours = neighbours/sum(weights(~isnan(neighbours)));
        neighbours((window^3-1)./2+1) = []; % Remove the center voxel
        if ~isempty(neighbours(~isnan(neighbours))) % Thus only excluding voxels with NaNs only as neighbors.
            NGTDM(value) = NGTDM(value) + abs(value-sum(neighbours(~isnan(neighbours))));
            countValid(value) = countValid(value) + 1;
        end
    end
end

end


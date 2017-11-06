function GLRLM = getGLRLM_M1(ROIonly,levels)

% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% - Xunkai Wei <xunkai.wei@gmail.com>
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
levelTemp = max(levels)+1;
ROIonly(isnan(ROIonly)) = levelTemp; % Last row needs to be taken out of the GLRLM
levels = [levels,levelTemp];


% QUANTIZATION EFFECTS CORRECTION
% In case (for example) we initially wanted to have 64 levels, but due to
% quantization, only 60 resulted.
uniqueVol = round(levels*adjust)/adjust;
ROIonly=round(ROIonly*adjust)/adjust;
NL = length(levels) - 1;


%INITIALIZATION
sizeV = size(ROIonly);
numInit = ceil(max(sizeV)*sqrt(3)); % Max run length
GLRLM = zeros(NL+1,numInit);


% START COMPUTATION
% Directions [1,0,0], [0 1 0], [1 1 0] and [-1 1 0] : 2D directions
% (x:right-left, y:top-bottom, z:3rd dimension)  
if numel(size(ROIonly)) == 3
    nComp = sizeV(3); % We can add-up the GLRLMs taken separately in every image in the x-y plane
else
    nComp = 1;
end
for i = 1:nComp
    image = ROIonly(:,:,i);
    uniqueIm = unique(image);
    NLtemp = length(uniqueIm);
    indexRow = zeros(NLtemp,1);
    temp = image;
    for j = 1:NLtemp
        indexRow(j) = find(uniqueIm(j)==uniqueVol);
        image(temp==uniqueIm(j)) = j;
    end
    
    % [1,0,0]
    GLRLMtemp = rle_0(image,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    
    % [0 1 0]
    GLRLMtemp = rle_0(image',NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    
    % [1 1 0]
    seq = zigzag(image);
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    
    % [-1 1 0]
    seq = zigzag(fliplr(image));
    GLRLMtemp = rle_45(seq,NLtemp);
    nRun = size(GLRLMtemp,2);
    GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
end

if numel(size(ROIonly)) == 3 % 3D DIRECTIONS
    % Directions [0,0,1], [1 0 1] and [-1 0 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp = sizeV(1); % We can add-up the GLRLMs taken separately in every image in the x-z plane
    image = zeros(sizeV(3),sizeV(2));
    for i = 1:nComp
        for j = 1:sizeV(3)
            image(j,1:end) = ROIonly(i,1:end,j);
        end
        uniqueIm = unique(image);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image;
        for j=1:NLtemp
            indexRow(j) = find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j)) = j;
        end
        
        % [0,0,1]
        GLRLMtemp = rle_0(image',NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        
        % [1 0 1]
        seq = zigzag(image);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        
        % [-1 0 1]
        seq = zigzag(fliplr(image));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    end

    % Directions [0,1,1] and [0 -1 1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    nComp = sizeV(2); % We can add-up the GLRLMs taken separately in every image in the y-z plane
    image = zeros(sizeV(1),sizeV(3));
    for i = 1:nComp
        for j = 1:sizeV(3)
            image(1:end,j) = ROIonly(1:end,i,j);
        end
        uniqueIm = unique(image);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image;
        for j = 1:NLtemp
            indexRow(j) = find(uniqueIm(j)==uniqueVol);
            image(temp==uniqueIm(j)) = j;
        end
        
        % [0,1,1]
        seq = zigzag(image);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        
        % [0 -1 1]
        seq = zigzag(fliplr(image));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    end

    % Four corners: [1,1,1], [-1,1,1], [-1,1,-1], [1,1,-1]
    % (x:right-left, y:top-bottom, z:3rd dimension)
    image = zeros(sizeV(3),sizeV(2));
    temp = rand(sizeV(3),sizeV(2));
    diagTemp = spdiags(temp);
    szDiag = size(diagTemp);
    diagMat1 = zeros(szDiag(1),szDiag(2),sizeV(1));
    diagMat2 = zeros(size(diagTemp,1),size(diagTemp,2),sizeV(1));
    for i = 1:sizeV(1)
        for j = 1:sizeV(3)
            image(j,1:end) = ROIonly(i,1:end,j);
        end
        try
            diagMat1(:,:,i)=spdiags(image);
        catch
            % Add a column at the beginning to prevent errors
            temp=spdiags(image);
            numberDiff=abs(size(temp,2)-size(diagMat1,2));
            if mod(numberDiff,2) % Odd difference number
                temp=padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat1(:,:,i)=temp(:,1:end-1);
            else
                diagMat1(:,:,i)=padarray(temp,[0,numberDiff/2,0],0);
            end
        end
        try
            diagMat2(:,:,i)=spdiags(fliplr(image));
        catch
            % Add a column at the beginning to prevent errors
            temp = spdiags(fliplr(image));
            numberDiff = abs(size(temp,2)-size(diagMat2,2));
            if mod(numberDiff,2) % Odd difference number
                temp = padarray(temp,[0,(numberDiff+1)/2,0],0);
                diagMat2(:,:,i) = temp(:,1:end-1);
            else
                diagMat2(:,:,i) = padarray(temp,[0,numberDiff/2,0],0);
            end
        end
    end
    for j = 1:szDiag(2)
        index = (diagMat1(:,j,1)~=0);
        nTemp = sum(index);
        image1 = zeros(sizeV(1),nTemp);
        image2 = zeros(sizeV(1),nTemp);
        for k = 1:sizeV(1)
            image1(k,1:nTemp) = diagMat1(index(1:end),j,k)';
            image2(k,1:nTemp) = diagMat1(index(1:end),j,k)';
        end
        
        % 2 first corners
        uniqueIm = unique(image1);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image1;
        for i = 1:NLtemp
            indexRow(i) = find(uniqueIm(i)==uniqueVol);
            image1(temp==uniqueIm(i)) = i;
        end
        seq = zigzag(image1);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        seq = zigzag(fliplr(image1));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        
        % 2 last corners
        uniqueIm = unique(image2);
        NLtemp = length(uniqueIm);
        indexRow = zeros(NLtemp,1);
        temp = image2;
        for i = 1:NLtemp
            indexRow(i) = find(uniqueIm(i)==uniqueVol);
            image2(temp==uniqueIm(i)) = i;
        end
        seq = zigzag(image2);
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
        seq = zigzag(fliplr(image2));
        GLRLMtemp = rle_45(seq,NLtemp);
        nRun = size(GLRLMtemp,2);
        GLRLM(indexRow(1:NLtemp),1:nRun) = GLRLM(indexRow(1:NLtemp),1:nRun) + GLRLMtemp(1:NLtemp,1:nRun); % Cumulative addition into the GLRLM
    end
end


% REMOVE UNECESSARY COLUMNS
GLRLM(end,:) = [];
stop = find(sum(GLRLM),1,'last');
GLRLM(:,(stop+1):end) = [];

end
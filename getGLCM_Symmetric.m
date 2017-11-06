function coocMat = getGLCM_Symmetric(varargin)
%inputStr = {TumorVolume,'Distance',[],'Direction',[],'numgray',levelsM+1};
%
%ljlubme@gmail.com
%Southern Medical University
%
%Default settings
coocMat= NaN;
distance = [1;2;4;8];
numLevels = 16;

offSet = [1 0 0; 1 1 0; 0 1 0; -1 1 0]; %2D Co-Occurrence directions 0,45,90,135 degrees
%the additional 9 directions of 3D volume
dimension3 = [0 0 1; 1 0 1; -1 0 1; 0,1,1; 0 -1 1; 1 1 1; -1 1 1; 1 1 -1; 1 1 -1];    
offSet = cat(1,offSet,dimension3);%13 directions

%checking inputs
data = varargin{1};
temp = size(data);
if size(temp)<3
    disp('Error: This program is designed for 3 dimensional data')
    return;
end
numInput = size(varargin,2);
for inputs =2:numInput
    temp = varargin{1,inputs};
    if ~ischar(temp)
        continue;
    end
    temp = upper(temp);
    switch (temp)
     
         case 'DIRECTION'
             temp2 = int8(varargin{1,inputs+1});
             if size(size(temp2),2) ~=2
                 disp('Error: Direction input is formatted poorly')
                 return;
             end
             if size(temp2,2) ~=3
                 disp(['Error: Incorrect number of columns in ' ... 
                     'direction variable'])
                 return;
             end
             if max(max(temp2))>1 | min(min(temp2))<-1
                 disp('Error: Direction values can only be {-1,0,1}')
                 return;
             end
             offSet = temp2;
        
        case 'DISTANCE'
            temp2 = int8(varargin{1,inputs+1});
            if size(size(temp2)) ~= 2
                disp('Error: Incorrect formatting of distance variable')
                return;
            end
            
            if sum(sum(size(temp2))) ~= max(size(temp2)+1)
                disp(['Error: Distance variable is to be a one ' ...
                    'dimensional array'])
                return;
            end
            
            distance = temp2;
            
        case 'NUMGRAY'
            temp2 = varargin{1,inputs+1};
            if temp2<1
                disp('The number of graylevels must be positive')
                return;
            end
            numLevels = uint16(temp2);
    end
end

noDirections = size(offSet,1); %number of directions, currently 13
coocMat = zeros(numLevels, numLevels, noDirections, size(distance,2));

for dist =1:size(distance,2) %distance
    [coocMat(:,:,:,dist)] = graycooc3d(data(:,:,:),distance(dist),numLevels,offSet);               
end
return

function [new_coMat]= graycooc3d(I,distance,numLevels,offSet)
%I = the 3D image matrix
%distance = a vector of the distances to analyze in
%numLevels = the number of graylevels to be used
%offSet = a matrix of the directions to analyze in
%coMat the Co-Occurrence matrices produced%%
%**************Variable initialization/Declaration**********************
%harMat =0;
noDirections = size(offSet,1); %number of directions, currently 13
coMat = zeros(numLevels,numLevels,noDirections);
%**************************Beginning analysis*************************
%Order of loops: Direction, slice, graylevel, graylevel locations
for direction =1:noDirections %currently 13 (for the 3d image)

    tempMat = zeros(numLevels,numLevels,size(I,3));
    for slicej =1:size(I,3)
         for j=1:numLevels %graylevel
             
             %find all the instances of that graylevel
            [rowj,colj] = find(I(:,:,slicej)==j);  

            %populating the Cooc matrix.
            for tempCount = 1:size(rowj,1) 
                rowT = rowj(tempCount) + distance*offSet(direction,1);
                colT = colj(tempCount) + distance*offSet(direction,2);
                sliceT = slicej + distance*offSet(direction,3);
                rowTnegative = rowj(tempCount) - distance*offSet(direction,1);%the symmetry of GLCM. 
                colTnegative = colj(tempCount) - distance*offSet(direction,2);%the symmetry of GLCM.
                sliceTnegative = slicej - distance*offSet(direction,3);%the symmetry of GLCM.      
                [I1, I2, I3] = size(I);  
                if rowT <= I1 && colT <= I2 && sliceT <= I3
                    if rowT > 0 && colT > 0 && sliceT > 0
                        %Error checking for NANs and Infinite numbers
                        IIntensity = I(rowT,colT,sliceT);
                       
                        if ~isnan(IIntensity)
                            if ~isinf(IIntensity)
                                tempMat(j,IIntensity,slicej)= tempMat...
                                    (j,IIntensity,slicej)+1;
                            end
                        end
                    end
                end
                  if rowTnegative <= I1 && colTnegative <= I2 && sliceTnegative <= I3
                    if rowTnegative > 0 && colTnegative > 0 && sliceTnegative > 0
                        
                        %Error checking for NANs and Infinite numbers
                        
                        IIntensitynegative = I(rowTnegative,colTnegative,sliceTnegative);% added by
                        if ~isnan(IIntensitynegative)
                            if  ~isinf(IIntensitynegative)
                                tempMat(j,IIntensitynegative,slicej)= tempMat...
                                    (j,IIntensitynegative,slicej)+1;
                                
                            end
                        end
                     end
                 end
            end
            
        end

    end
    for slicej =1:size(I,3)
        coMat(:,:,direction)= coMat(:,:,direction)+tempMat(:,:,slicej);
    end   
end
new_coMat=coMat(:,:,:);
return




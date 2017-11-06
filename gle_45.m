function oneglglm = gle_45(seq,NL)
% RLE   image gray level Gap Length matrix for 45 and 135
% This file is to handle the zigzag scanned sequence for 45 or 135 degree
% direction. Note for 135, just swap the left and the right colum
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
%  -------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007


% Assure row number is exactly the gray leve;
% number of seqence
m =length(seq);
% number to store the possible max coloums
n = findmaxnum(seq);
%

oneglglm=zeros(NL,n);

for i=1:m
    x=seq{i};
    % run length Encode of each vector
    for j = 1:max(x(:))
        index = find(x == j);
        if ~isempty(index) && length(index)~=1
            for k = 1:length(index)-1
                for l = k+1:length(index)
                gap = index(l)-index(k);
                oneglglm(j,gap) = oneglglm(j,gap) +1;
                end
            end
        else
        end
    end   
end
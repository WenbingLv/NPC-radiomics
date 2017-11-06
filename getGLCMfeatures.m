function features = getGLCMfeatures(GLCM)
%
%ljlubme@gmail.com
%Southern Medical University
%
% PRELIMINARY
features = struct;
matrixtemp = GLCM;
GLCM = GLCM/(sum(GLCM(:))); % Normalization of GLCM
nL = max(size(GLCM));
indVect = 1:nL;
[colGrid,rowGrid] = meshgrid(indVect,indVect);


% COMPUTATION OF TEXTURE FEATURES
% 1. Energy
features.Energy = sum(sum(GLCM.^2));

% 2. Entropy
features.Entropy = -sum(sum(GLCM.*log(GLCM + realmin)));

% 3. Difference entropy
temp = zeros(1,nL);
for k = 0:nL-1
  for i = 1:nL
     for j = 1:nL
       if(abs(i-j)==k)
        temp(k+1) = temp(k+1) + GLCM(i,j);
        else
        end
     end
  end
end
px_y = temp;
features.DiffEntropy = -px_y*(log(px_y+realmin))';

%4. Sum entropy
temp = zeros(1,2*nL-2+1);
for k = 2:2*nL
  for i = 1:nL
     for j = 1:nL
       if(i+j==k)
        temp(k-1) = temp(k-1) + GLCM(i,j);
        else
        end
     end
  end
end
pxy = temp;
features.SumEntropy = -pxy*(log(pxy+realmin))';

% 5. Variance and 24. SumAverage1
ui = indVect*sum(GLCM,2);
uj = indVect*sum(GLCM)';
tempS = rowGrid.*GLCM + colGrid.*GLCM;
tempV = (rowGrid-ui).^2.*GLCM + (colGrid-uj).^2.*GLCM;
features.SumAverage1 = sum(tempS(:))/(nL^2);% modified
features.Variance = sum(tempV(:))/(nL^2);%equal: 0.5 is unnecessary
% features.SumAverage = 0.5*sum(tempS(:))/(nL^2);% origional code,0.5
% features.Variance = 0.5*sum(tempV(:))/(nL^2);

% 6. Sum of squares variance
u = mean(GLCM(:));
temp = 0;
for i = 1:nL;
    for j = 1:nL;
        temp = temp + (i-u)^2*GLCM(i,j);
    end
end
features.SumSquVar = temp;

% 8. Max possibility
features.MaxPossibility = max(max(GLCM));

% 9. Contrast
contrast = 0.0;
for n = 0:nL-1%sum start from 1 based eq contrast
   temp = 0;
   for i = 1:nL
      for j = 1:nL
         if (abs(i-j) == n)
            temp = temp+GLCM(i,j);
         end
      end
   end
   contrast = contrast + n^2*temp;
end
features.Contrast = contrast;

% 10. Dissimilarity
diffMat = abs(rowGrid-colGrid);
temp = diffMat.*GLCM;
features.Dissimilarity = sum(temp(:));

% 11. Homogeneity
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + GLCM(i,j)/(1+abs(i-j));
   end
end
features.Homogeneity = temp;

% 12. Inverse Different Moment
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + GLCM(i,j)/(1+abs(i-j).^2);
   end
end
features.InDiffMoment = temp;% also called local homogeneity

% 13. Correlation
features.Correlation = graycoprops(round(matrixtemp),'Correlation');
features.Correlation = struct2cell(features.Correlation);
features.Correlation = features.Correlation{1};

% 14. Difference Variance
ux_y = indVect*px_y';
diffvar = -(indVect-ux_y).^2*px_y';
features.DiffVar = diffvar;

%  15. Autocorrelation
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + i*j*GLCM(i,j);
   end
end
features.AutoCorrelation = temp;

%  16. Cluster prominence
ux = indVect*sum(GLCM,2);
uy = indVect*sum(GLCM)';
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux-uy).^4.*GLCM(i,j);
   end
end
features.ClusterPro = temp;

% 17. CLuster shade
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux-uy).^3.*GLCM(i,j);
   end
end
features.ClusterShade = temp;

% 18. Cluster tendency
temp = 0;
for i = 1:nL
   for j = 1:nL
      temp = temp + (i+j-ux-uy).^2.*GLCM(i,j);
   end
end
features.ClusterTen = temp;

% 19. 20. informational measure of correlation1(ICM1and ICM2)
H = features.Entropy;
px = sum(GLCM,2);%sum of the row
py = sum(GLCM);%sum of the clumn
HXY1 = 0;
HXY2 = 0;
for i = 1:nL
   for j = 1:nL
      HXY1 = HXY1-GLCM(i,j)*log(px(i)*py(j)+realmin);
      HXY2 = HXY2-px(i)*py(j)*log(px(i)*py(j)+realmin);
   end
end
HX = -px'*(log(px+realmin));%entropies of px and py
HY = -py*(log(py+realmin))';
features.ICM1 = (H-HXY1)./(max(HX,HY)+realmin);
features.ICM2 = (1-exp(-2*(HXY2-H))).^0.5;

% 21. Inverse variance 22. Inverse difference moment normalized
tempIva = 0;
tempIDMN = 0;
tempIDN = 0;
for i = 1:nL
   for j = 1:nL
       tempIDMN = tempIDMN + GLCM(i,j)/(1+(i-j)^2/(nL^2));
       tempIDN = tempIDN + GLCM(i,j)/(1+(i-j)/nL);
       if i~=j
          tempIva = tempIva + GLCM(i,j)/((i-j)^2);
       else
       end
   end
end
features.InVar = tempIva;% actually it is Inverse variance,should  change to InVar
features.IDMN = tempIDMN;

% 23. IDN
features.IDN = tempIDN;

% 25. Sum average2
i = 2:2*nL;
features.SumAverage2 = i*pxy';

%7. Sum variance
SA = features.SumAverage2;
temp = 0;
for i = 1:length(pxy)
    temp = temp + (i+1-SA).^2*pxy(i);
end
features.SumVar = temp;

% 26. Agreement
pe = 0;
p0 = trace(GLCM);
Ag = 0;
for i = 1:nL;
    pe = pe + GLCM(i,:)*GLCM(:,i);
end
Ag = (p0 - pe)/(1-pe+realmin);
features.Agreement = Ag;

end









function features = getGLSZMfeatures(GLSZM)
% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% - Revision: Oct 2017
% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLSZM); % Size of GLSZM
nRuns = sum(GLSZM(:));
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLSZM
pg = sum(GLSZM,2)'; 
pr = sum(GLSZM); 

% COMPUTATION OF TEXTURE FEATURES
% 1. Small Zone Emphasis (SZE)
features.SZE = (pr*(cVect.^(-2))')/nRuns;

% 2. Large Zone Emphasis (LZE)
features.LZE = (pr*(cVect.^2)')/nRuns;

% 3. Gray-Level Nonuniformity (GLN)
features.GLN = sum(pg.^2)/nRuns;

% 4. Zone-Size Nonuniformity (ZSN)
features.ZSN = sum(pr.^2)/nRuns;

% 5. Zone Percentage (ZP)
features.ZP = nRuns/(pr*cVect');

% 6. Low Gray-Level Zone Emphasis (LGZE)
features.LGZE = (pg*(rVect.^(-2))')/nRuns;

% 7. High Gray-Level Zone Emphasis (HGZE)
features.HGZE = (pg*(rVect.^2)')/nRuns;

% 8. Small Zone Low Gray-Level Emphasis (SZLGE)
features.SZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^(-2))))/nRuns;

% 9. Small Zone High Gray-Level Emphasis (SZHGE)
features.SZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^(-2))))/nRuns;

% 10. Large Zone Low Gray-Level Emphasis (LZLGE)
features.LZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^2)))/nRuns;

% 11. Large Zone High Gray-Level Emphasis (LZHGE)
features.LZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^2)))/nRuns;


% New features
pg = sum(GLSZM,2)'; pr = sum(GLSZM);
ug = (pg*rVect')/(sz(1)*sz(2));
ur = (pr*cVect')/(sz(1)*sz(2));

% 12. Gray-Level Variance (GLV)
GLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GLV = GLV + (GLSZM(g,r)*g-ug)^2;
    end
end
features.GLV = GLV/(sz(1)*sz(2));

% 13. Zone-Size Variance (ZSV)
ZSV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        ZSV = ZSV + (GLSZM(g,r)*r-ur)^2;
    end
end
features.ZSV = ZSV/(sz(1)*sz(2));

end
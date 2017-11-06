function features = getGLRLMfeatures(GLRLM)
% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% - Revision: Oct 2017

% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLRLM); % Size of GLRLM
nRuns = sum(GLRLM(:));
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLRLM
pg = sum(GLRLM,2)'; % Gray-Level Run-Number Vector
pr = sum(GLRLM); % Run-Length Run-Number Vector


% COMPUTATION OF TEXTURE FEATURES
% 1. Short Run Emphasis (SRE)
features.SRE = (pr*(cVect.^(-2))')/nRuns;

% 2. Long Run Emphasis (LRE)
features.LRE = (pr*(cVect.^2)')/nRuns;

% 3. Gray-Level Nonuniformity (GLN)
features.GLN = sum(pg.^2)/nRuns;

% 4. Run-Length Nonuniformity (RLN)
features.RLN = sum(pr.^2)/nRuns;

% 5. Run Percentage (RP)
features.RP = nRuns/(pr*cVect');

% 6. Low Gray-Level Run Emphasis (LGRE)
features.LGRE = (pg*(rVect.^(-2))')/nRuns;

% 7. High Gray-Level Run Emphasis (HGRE)
features.HGRE = (pg*(rVect.^2)')/nRuns;

% 8. Short Run Low Gray-Level Emphasis (SRLGE)
features.SRLGE = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^(-2))))/nRuns;

% 9. Short Run High Gray-Level Emphasis (SRHGE)
features.SRHGE = sum(sum(GLRLM.*(rMat.^2).*(cMat.^(-2))))/nRuns;

% 10. Long Run Low Gray-Level Emphasis (LRLGE)
features.LRLGE = sum(sum(GLRLM.*(rMat.^(-2)).*(cMat.^2)))/nRuns;

% 11. Long Run High Gray-Level Emphasis (LRHGE)
features.LRHGE = sum(sum(GLRLM.*(rMat.^2).*(cMat.^2)))/nRuns;


% New features
GLRLM = GLRLM./nRuns;
pg = sum(GLRLM,2)'; pr = sum(GLRLM);
ug = (pg*rVect')/(sz(1)*sz(2));
ur = (pr*cVect')/(sz(1)*sz(2));

% 12. Gray-Level Variance (GLV)
GLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GLV = GLV + (GLRLM(g,r)*g-ug)^2;
    end
end
GLV = GLV/(sz(1)*sz(2));
features.GLV = GLV;

% 13. Run-Length Variance (RLV)
RLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        RLV = RLV + (GLRLM(g,r)*r-ur)^2;
    end
end
RLV = RLV/(sz(1)*sz(2));
features.RLV = RLV;

end
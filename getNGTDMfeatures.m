function features = getNGTDMfeatures(NGTDM,countValid)
% AUTHOR(S): 
% - Lijun Lu <ljlubme@gmail.com>
% - Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% - Revision: Oct 2017
% PRELIMINARY
nTot = sum(countValid);
countValid = countValid./nTot; 
NL = length(NGTDM);
Ng = sum(countValid~=0);
pValid = find(countValid>0);
nValid = length(pValid);


% COMPUTATION OF TEXTURES
% 1. Coarseness
features.Coarseness = (((countValid')*NGTDM) + eps)^(-1);

% 2. Contrast
val = 0;
for i = 1:NL
    for j = 1:NL
        val = val + countValid(i)*countValid(j)*(i-j)^2;
    end
end
features.Contrast = val*sum(NGTDM)/(Ng*(Ng-1)*nTot);

% 3. Busyness
denom = 0;
for i = 1:nValid
    for j = 1:nValid
        denom = denom + abs(pValid(i)*countValid(pValid(i))-pValid(j)*countValid(pValid(j)));
    end
end
features.Busyness = ((countValid')*NGTDM)/denom;

% 4. Complexity
val = 0;
for i = 1:nValid
    for j = 1:nValid
        val = val + (abs(pValid(i)-pValid(j))/(nTot*(countValid(pValid(i)) + countValid(pValid(j)))))*(countValid(pValid(i))*NGTDM(i) + countValid(pValid(j))*NGTDM(j));
    end
end
features.Complexity = val;

% 5. Strength
val = 0;
for i = 1:nValid
    for j = 1:nValid
        val = val + (countValid(pValid(i))+countValid(pValid(j)))*(pValid(i)-pValid(j))^2;
    end
end
features.Strength = val/(eps+sum(NGTDM));

end
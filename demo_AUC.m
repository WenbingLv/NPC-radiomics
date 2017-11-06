%this is an example for investigating the classification performance of LGRE extracted by
%strategy M1 and M13, logistic regression and leave-one-out corss validation are used,
%the AUC, sensitivity and specificity are also calculated.
%
%ljlubme@gmail.com
%Southern Medical University
%
clear;clc;

load LGRE_M1_M13.mat
load ALL_label.mat
Feature = LGRE_M1_M13;% Low Gray Level Run Emphasis (LGRE) feature, 106 cases * 2 extraction methods (M1 and M13)
id = 1:length(ALL_label);% 69 cases of NPC with label of 1, 37 cases of CN with label of 0
for i = 1:2%i=1-->M1, i=2-->M13, 
      for j = 1:length(ALL_label)%j, the number of cases
             idTrain = setdiff(id,j);%leave-one-out corss validation
             b = glmfit(Feature(idTrain,i),ALL_label(idTrain),'binomial');  % logistic regression
             p(j,:) = glmval(b,Feature(j,i),'logit');%out put probability
     end
    [X,Y,T,AUCa] = perfcurve(ALL_label,p,1);% X and Y indicate fpr and tpr respevtively, roc is equal to perfcurve            
    AUC_all(i) = AUCa;
    Z = Y+(1-X);% find out the point in ROC that maximized the sum of sensitvity and specificity
    idok = find(Z==max(Z));
    Sen_all(i) = Y(idok);%sensitivity
    Spe_all(i) = 1-X(idok);%specificity
    T_all(i) = T(idok);%threshold
    clear p;
end 

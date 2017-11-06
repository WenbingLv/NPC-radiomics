%this is an example of how to calculate the ICC of LGRE extracted by
%strategy M1 and M13.
%
%ljlubme@gmail.com
%Southern Medical University
%
clear;clc;

load LGRE_M1_M13.mat
Feature = LGRE_M1_M13;% Low Gray Level Run Emphasis (LGRE) feature, 106 cases * 2 extraction methods (M1 and M13)
[p,anovatab,stats] = anova1(Feature,[],'off'); %one-way ANOVA
BMS = anovatab{2,4};%between_subject mean squares
WMS = anovatab{3,4};%within_subject mean squares
ICC = (BMS - WMS)./(BMS + WMS); 

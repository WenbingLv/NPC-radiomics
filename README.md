# NPC-radiomics
This is a matlab demo showing how to construct radiomic matrix and extract radiomic features using different strategies/parameter 
settings for PET/CT scanned NPC images. 

Citation: 

[1] Lv W, Yuan Q et al. Impact of parameter settings for radiomic features on their robustness and diagnostic task performance: application to nasopharyngeal 18F-FDG PET/CT. European Radiology, under review.

[2] Lu L, Lv W, Jiang J et al (2016) Robustness of radiomic features in [11C]choline and [18F]FDG PET/CT imaging of nasopharyngeal carcinoma: impact of segmentation and discretization. Mol Imaging Biol 18:935-945

Author: Wenbing Lv, Lijun Lu  <ljlubme@gmail.com>
Southern Medical University

History:

Creation: 2017.11

This matlab code contains 3 deoms:

1. “demo_FeaExtr”: matlab codes to load PET images, segmentation mask and extract radiomic features under different parameter settings, i.e. symmetry, averaging strategy and distance for GLCM,  averaging strategy for GLRLM, neighborhood extent for GLSZM and window size for NGTDM.

2. "demo_ICC":matlab codes to calculate the intra-class correlation coeffcient of feature extracted under different parameter settings.

3. "demo_AUC":matlab codes to calculate the AUC, Sensitivity and Specificity after logistic regression with leave-one-out cross validation.

Ackonwledgements:
- Martin Vallières: https://github.com/mvallieres/radiomics
- Wei's GLRLM toolbox: Xunkai Wei, Gray Level Run Length Matrix Toolbox
  v1.0, Software,Beijing Aeronautical Technology Research Center, 2007.
  <http://www.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox>
- Q. Li: <http://www.mathworks.com/matlabcentral/fileexchange/23377-ellipsoid-fitting>
- CERR development team: <http://www.cerr.info/>
- Dirk-Jan Kroon (imresize3D.m): <http://www.mathworks.com/matlabcentral/fileexchange/21451-multimodality-non-rigid-demon-algorithm-image-registration/content//functions/imresize3d.m>
- David Reshef and Yakir Reshef: MINE version 1.0.1d <http://www.exploredata.net/> 
- DREES development team: <http://www.cerr.info/drees>
- Enric Junqué de Fortuny (fastAUC.cpp): <http://www.mathworks.com/matlabcentral/fileexchange/41258-faster-roc-auc>
- François Beauducel (roundsd.m): <http://www.mathworks.com/matlabcentral/fileexchange/26212-round-with-significant-digits>
- Jos van der Geest (herrorbar.m): <http://www.mathworks.com/matlabcentral/fileexchange/3963-herrorbar> 

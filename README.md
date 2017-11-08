# NPC-radiomics
This is a matlab demo showing how to construct radiomic matrix and extract radiomic features using different strategies/parameter 
settings for PET/CT scanned NPC images. 

Citation: 

[1] W.Lv, Q.Yuan et al. Impact of parameter settings for radiomic features on their robustness and diagnostic task performance: application to nasopharyngeal 18F-FDG PET/CT. European Radiology,under review.

[2] Lu L, Lv W, Jiang J et al (2016) Robustness of radiomic features in [11C]choline and [18F]FDG PET/CT imaging of nasopharyngeal carcinoma: impact of segmentation and discretization. Mol Imaging Biol 18:935-945

Reference:

[1] Vallières, M. et al. (2015). A radiomics model from joint FDG-PET and MRI texture features for the prediction of lung metastases in soft-tissue sarcomas of the extremities. Physics in Medicine and Biology, 60(14), 5471-5496. doi:10.1088/0031-9155/60/14/5471

[2] Zhou, H., Vallieres, M., Bai, H.X. et al. (2017). MRI features predict survival and molecular markers in diffuse lower-grade gliomas. Neuro-Oncology, XX(XX), 1-10. doi:10.1093/neuonc/now256

Author: Wenbing Lv, Lijun Lu <ljlubme@gmail.com>
Southern Medical University

History:
- Creation: 2017.11

This matlab code antains 3 deoms:

1. “demo_FeaExtr”: matlab codes to load PET images, segmentation mask and extract radiomic features under different parameter settings, i.e. symmetry, averaging strategy and distance for GLCM,  averaging strategy for GLRLM, neighborhoods extent for GLSZM and window size for NGTDM.

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


            using two pairs of fine and coarse images
;         the program can be used for whole tm scene and vi index product
;   a fast version of estarfm code: improve the effeciency of original estarfm 
;   through controling the number of similar pixels
;        developed by (1) zhu xiaolin,email: zhuxiaolin55@gmail.com
;             department of land surveying and geo-informatics
;             the hong kong polytechnic university
;             
;            debugging history: 
;            1)5/12/2012 correct the abnormal prediction
;            2)9/29/2013 correct the spatial distance calculation for integrate window 
;            3)7/10/2014 correct the abnormal value of spectral distance and use all bands to indentify background
;            4)2/13/2017 add one parameter to specify the value of pixels of background or missing
;            5)1/1/2018  improve efficiency and modified the weight for fusing vi index
;            6)3/11/2008 correct a bug in spatial distance caculation
;            7)7/27/2018 correct a bug of abnormal coversion coefficent estimation when the two input pairs are too similar
;            8)7/27/2018 improve the prediction when no enough simular pixels are selected
;                        
;please cite the reference: xiaolin zhu, jin chen, feng gao, & jeffrey g masek.
;an enhanced spatial and temporal adaptive reflectance fusion model for complex
;heterogeneous regions. remote sensing of environment,2010,114,2610-2623
;
;                     copyright belongs to xiaolin zhu


YOU can download from website of https://xiaolinzhu.weebly.com/open-source-code.html

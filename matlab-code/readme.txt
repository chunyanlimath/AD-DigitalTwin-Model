This folder includes a template MATLAB code for completing the hierarchical structure and homotopy regularization technique for the Abeta biomarker. And all of the other biomarkers can be computed in a similar manner.  
Author: Chunyan Li
Update: 5/13/25

1. Population mean of FC matrix is stored in 'Lp68.mat' and the population mean of Adjacent matrix is saved in 'A.mat'
  2. homogeneous model parameters and i.c. 
    3. collect scalar results
      4. 68D model parameters 
         5. collect 68D results
           6. inference u0 (i.c.) 
             7. collect i.c. results
               8. learn L = Lp + f(u,v)
                 9. collect results.
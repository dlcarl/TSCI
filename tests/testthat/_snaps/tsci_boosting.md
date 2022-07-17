# snapshot tsci_boosting

    Code
      output
    Output
      
      Statistics about the data splitting procedure:
      Sample size A1: 67 
      Sample size A2: 33 
      Number of data splits: 2 
      Aggregation method: FWER 
      
      Coefficients:
                  Estimate Std_Error   2.5 %  97.5 % Pr(>|t|)
      TSCI-comp    1.10005         . 0.93855 1.26433        0
      TSCI-robust  1.10005         . 0.93855 1.26433        0
      OLS          1.27239         . 1.15407 1.39072        0
      TSCI-Cor-q0  1.20214         . 1.03944 1.36782        0
      TSCI-Cor-q1  1.10005         . 0.93855 1.26433        0
      TSCI-Cor-q2  1.20110         . 0.95743 1.43599        0
      TSCI-Cor-q3  1.04459         . 0.75894 1.36721        0
      TSCI-Cor-q4  1.20654         . 0.84331 1.51674        0
      
      Statistics about the treatment model:
      Estimation method: L2 Gradient Tree Boosting 
      Residual standard error: 1.8963 
      R-squared: 0.5925 
      
      Statistics about the outcome model:
                  Residual_Standard_Error R_Squared
      OLS                          0.9503    0.8557
      TSCI-Cor-q0                  0.9408    0.8584
      TSCI-Cor-q1                  0.8834    0.8770
      TSCI-Cor-q2                  0.8389    0.8908
      TSCI-Cor-q3                  0.8670    0.8852
      TSCI-Cor-q4                  0.8332    0.8956
      
      Statistics about the violation space selection:
          q_comp q_robust Qmax
      OLS      0        0    0
      q0       0        0    0
      q1       2        2    2
      q2       0        0    0
      q3       0        0    0
      q4       0        0    0
      
      Statistics about the IV strength:
         IV_Strength IV_Threshold
      q0       95.06        65.74
      q1       84.73        58.55
      q2       23.50        46.22
      q3       14.84        44.03
      q4       10.78        44.98


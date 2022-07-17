# snapshot tsci_forest

    Code
      output
    Output
      
      Statistics about the data splitting procedure:
      Sample size A1: 67 
      Sample size A2: 33 
      Number of data splits: 2 
      Aggregation method: FWER 
      
      Coefficients:
                  Estimate Std_Error    2.5 %  97.5 % Pr(>|t|)
      TSCI-comp    1.04587         .  0.82594 1.25453  0.00000
      TSCI-robust  1.04587         .  0.82594 1.25453  0.00000
      OLS          1.27239         .  1.15407 1.39072  0.00000
      TSCI-Cor-q0  1.19371         .  1.00566 1.38104  0.00000
      TSCI-Cor-q1  1.04587         .  0.82594 1.25453  0.00000
      TSCI-Cor-q2  1.00225         .  0.27357 1.78508  0.00369
      TSCI-Cor-q3  0.75068         . -0.37493 1.80200  0.26916
      TSCI-Cor-q4  0.85212         . -1.14394 2.16635  0.85597
      
      Statistics about the treatment model:
      Estimation method: Random Forest 
      Residual standard error: 1.9154 
      R-squared: 0.4975 
      
      Statistics about the outcome model:
                  Residual_Standard_Error R_Squared
      OLS                          0.9503    0.8557
      TSCI-Cor-q0                  1.0006    0.8530
      TSCI-Cor-q1                  0.9582    0.8670
      TSCI-Cor-q2                  0.9191    0.8791
      TSCI-Cor-q3                  1.0655    0.8375
      TSCI-Cor-q4                  1.1215    0.8143
      
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
      q0       87.36        46.03
      q1       72.47        37.29
      q2        1.69        14.39
      q3        1.23        13.23
      q4        0.28        10.72


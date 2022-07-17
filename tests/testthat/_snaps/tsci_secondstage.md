# snapshot tsci_secondstage

    Code
      output
    Output
      
      Statistics about the data splitting procedure:
      Sample size: 100 
      No sample splitting was performed.
      
      Coefficients:
                  Estimate Std_Error    2.5 %   97.5 % Pr(>|t|)
      TSCI-comp    1.06057   0.06921  0.92493  1.19622  0.00000
      TSCI-robust  1.06057   0.06921  0.92493  1.19622  0.00000
      OLS          1.27239   0.05279  1.16893  1.37586  0.00000
      TSCI-Cor-q0  1.15941   0.06846  1.02524  1.29358  0.00000
      TSCI-Cor-q1  1.06057   0.06921  0.92493  1.19622  0.00000
      TSCI-Cor-q2  9.95152   3.64748  2.80260 17.10045  0.72285
      TSCI-Cor-q3  0.66071   3.54733 -6.29193  7.61335  0.79429
      TSCI-Cor-q4  1.33479   0.23375  0.87665  1.79293  0.00000
      
      Statistics about the treatment model:
      Estimation method: Specified by User 
      Residual standard error: 1.771 
      R-squared: 0.5089 
      
      Statistics about the outcome model:
                  Residual_Standard_Error R_Squared
      OLS                          0.9503    0.8557
      TSCI-Cor-q0                  0.9644    0.8498
      TSCI-Cor-q1                  0.9096    0.8678
      TSCI-Cor-q2                  2.4195    0.0739
      TSCI-Cor-q3                  1.0381    0.8313
      TSCI-Cor-q4                  0.8093    0.8985
      
      Statistics about the violation space selection:
          q_comp q_robust Qmax
      OLS      0        0    0
      q0       0        0    0
      q1       1        1    1
      q2       0        0    0
      q3       0        0    0
      q4       0        0    0
      
      Statistics about the IV strength:
         IV_Strength IV_Threshold
      q0      184.57        78.36
      q1      160.74        72.88
      q2        0.50        22.67
      q3        0.00        10.00
      q4        0.00        10.00


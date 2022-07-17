# snapshot tsci_poly

    Code
      output
    Output
      
      Statistics about the data splitting procedure:
      Sample size: 100 
      No sample splitting was performed.
      
      Coefficients:
                  Estimate Std_Error   2.5 %  97.5 % Pr(>|t|)
      TSCI-comp    1.05508   0.06770 0.92239 1.18776  0.00000
      TSCI-robust  1.05508   0.06770 0.92239 1.18776  0.00000
      OLS          1.27239   0.05279 1.16893 1.37586  0.00000
      TSCI-Cor-q0  1.16431   0.06633 1.03432 1.29431  0.00000
      TSCI-Cor-q1  1.05508   0.06770 0.92239 1.18776  0.00000
      TSCI-Cor-q2  1.14193   0.13072 0.88573 1.39814  0.00000
      TSCI-Cor-q3  0.76652   0.30467 0.16938 1.36366  0.00247
      TSCI-Cor-q4  1.13063   0.46217 0.22479 2.03647  0.01017
      TSCI-Cor-q5  5.10575   1.49000 2.18540 8.02610  0.55221
      
      Statistics about the treatment model:
      Estimation method: OLS with Polynomials 
      Residual standard error: 1.7727 
      R-squared: 0.5182 
      
      Statistics about the outcome model:
                  Residual_Standard_Error R_Squared
      OLS                          0.9503    0.8557
      TSCI-Cor-q0                  0.9615    0.8507
      TSCI-Cor-q1                  0.9098    0.8677
      TSCI-Cor-q2                  0.8628    0.8822
      TSCI-Cor-q3                  1.0006    0.8432
      TSCI-Cor-q4                  0.8341    0.8922
      TSCI-Cor-q5                  2.0585    0.3504
      
      Statistics about the violation space selection:
          q_comp q_robust Qmax
      OLS      0        0    0
      q0       0        0    0
      q1       1        1    1
      q2       0        0    0
      q3       0        0    0
      q4       0        0    0
      q5       0        0    0
      
      Statistics about the IV strength:
         IV_Strength IV_Threshold
      q0      195.26        78.02
      q1      172.17        72.63
      q2       31.21        32.16
      q3       13.14        34.39
      q4        3.29        22.88
      q5        0.00        10.00


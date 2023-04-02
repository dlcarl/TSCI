# TSCI 2.0.0

## breaking changes

* Bootstrap options for treatment effect standard error estimation and
  violation space selection were added to `tsci_boosting`, `tsci_forest`, 
  `tsci_poly` and `tsci_secondstage`. These bootstrap options are set to be the default method.
  
* An upper limit of 40 for the threshold for the IV strength was added to `tsci_boosting`, `tsci_forest`, 
  `tsci_poly` and `tsci_secondstage`.
  
## bug fixes

* Fixed an issue that led to presenting wrong p-values when using `tsci_poly` and `tsci_secondstage`.

* Fixed an issue that led to an error when trying to use parallelization on Windows.


# TSCI 1.0.1

* First version of the package. 

# TSCI 3.0.4

## Minor improvements and fixes

* Performance optimization

# TSCI 3.0.3

## Minor improvements and fixes

* Limited the random forest and boosting algorithms to only use one thread.
* Modified examples and tests to decrease computational burden.


# TSCI 3.0.2

## Minor improvements and fixes

* Modified examples and tests to decrease computational burden.

# TSCI 3.0.1

## Minor improvements and fixes

* Added information about reproducibility issue of `tsci_forest()' when using     different operating systems to the documentation.

# TSCI 3.0.0

## breaking changes

* `coef.tsci()` does now accept the argument `parm`, replacing `all`, that will provide the same
  functionality as in `confint.tsci()`. 
  
## Minor improvements and fixes

* `summary.tsci()` will now print p-values as by `format.pval()`.

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

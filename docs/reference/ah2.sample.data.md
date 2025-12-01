# Generate a sample data from the pbc data

This is a function to retrieve 312 randomized patients from the pbc data
in survival package.

## Usage

``` r
ah2.sample.data()
```

## Details

The function creates a sample dataset to illustrate the usage of the
function [`ah2()`](https://www.uno1lab.com/survAH/reference/ah2.md) in
this package. The original pbc data in `survival` package consists of
418 patients data. This function loads the pbc data, select the 312
patients who were randomized. The status variable is edited, so that 1
indicates death and 0 indicates alive.

## See also

`pbc` in survival package

## Examples

``` r
D=surv2test.sample.data()
#> Error in surv2test.sample.data(): could not find function "surv2test.sample.data"
head(D)
#>                               
#> 1 function (expr, name)       
#> 2 .External(C_doD, expr, name)
```

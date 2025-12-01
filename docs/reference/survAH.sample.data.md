# Generate a sample data from the pbc data

This is a function to retrieve 312 randomized patients from the pbc data
in survival package.

## Usage

``` r
survAH.sample.data()
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
D=survAH.sample.data()
head(D)
#>        time status arm
#> 1  1.095140      1   1
#> 2 12.320329      0   1
#> 3  2.770705      1   1
#> 4  5.270363      1   1
#> 5  4.117728      0   0
#> 6  6.852841      1   0
```

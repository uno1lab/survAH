# Generate a sample data from the pbc data

This is a function to retrieve 312 randomized patients from the pbc
data.

## Usage

``` r
ahreg.sample.data(t.unit="year")
```

## Arguments

- t.unit:

  Specify the time unit. It supports "year" (default), "month", and
  "day".

## Value

returns a data frame

## Details

The function creates a sample dataset to illustrate the usage of the
function [`ahreg()`](https://www.uno1lab.com/survAH/reference/ahreg.md)
in this package. The original pbc data in `survival` package consists of
418 patients data. This function loads the pbc data, select the 312
patients who were randomized. The status variable is edited, so that 1
indicates death and 0 indicates alive.

## See also

`pbc` in survival package

## Examples

``` r
D=ahreg.sample.data()
head(D)
#>        time status arm      age edema bili albumin protime
#> 1  1.095140      1   1 58.76523   1.0 14.5    2.60    12.2
#> 2 12.320329      0   1 56.44627   0.0  1.1    4.14    10.6
#> 3  2.770705      1   1 70.07255   0.5  1.4    3.48    12.0
#> 4  5.270363      1   1 54.74059   0.5  1.8    2.54    10.3
#> 5  4.117728      0   0 38.10541   0.0  3.4    3.53    10.9
#> 6  6.852841      1   0 66.25873   0.0  0.8    3.98    11.0
```

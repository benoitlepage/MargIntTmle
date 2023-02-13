
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MargIntTmle

<!-- badges: start -->
<!-- badges: end -->

The goal of MargIntTmle is to estimate marginal interaction effects
using g-computation, IPTW or TMLE.

## Installation

You can install the development version of MargIntTmle from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("benoitlepage/MargIntTmle")
```

## First example

In this first example, we simulate a data set of `N` = 1000 rows, with
three baseline confounders (`conf1`, `conf2`, and `conf3`), two
exposures (`sex` and `env`) and one outcome `hlth.outcome`. We used the
default parameters defined in the `param.causal.model()` function.

``` r
require(MargIntTmle)
#> Le chargement a nécessité le package : MargIntTmle

set.seed(12345)
df <- generate.data(N = 1000, b = param.causal.model())
head(df)
#>   conf1 conf2 conf3 sex env hlth.outcome
#> 1     1     0     0   0   0            0
#> 2     1     1     1   1   0            0
#> 3     1     0     1   0   0            0
#> 4     1     0     0   1   0            0
#> 5     0     0     1   0   0            0
#> 6     0     0     0   0   0            0
```

The `int.ltmleMSM` function is used to call the `ltmleMSM` function from
the `ltmle` package, in order to estimate marginal interaction effects
from a Marginal structural model. In this example, we use the TMLE
estimator.

Several quantities of interest are than calculated using the
`estim.int.effects` function.

``` r
require(ltmle)
require(SuperLearner)

# define Q and g formulas following the argument notation of the ltmle package:
Q_formulas = c(hlth.outcome="Q.kplus1 ~ conf1 + conf2 + conf3 + sex * env")
g_formulas = c("sex ~ conf1 + conf2","env ~ conf1 + conf3")

# choose a set of fitting libraries to pass to SuperLearner:
SL.library = list(Q=list("SL.glm"),g=list("SL.glm"))

# apply the int.ltmleMSM function. In order to apply the TMLE and IPTW estimators, 
# gcomp argument is set to FALSE.
interaction.ltmle <- int.ltmleMSM(data = df,
                                  Qform = Q_formulas,
                                  gform = g_formulas,
                                  Anodes = c("sex", "env"),
                                  Lnodes = c("conf1", "conf2", "conf3"),
                                  Ynodes = c("hlth.outcome"),
                                  SL.library = SL.library,
                                  gcomp = FALSE,
                                  iptw.only = FALSE,
                                  survivalOutcome = FALSE,
                                  variance.method = "ic")

# several quantities of interest for interaction effects are calculated using the 
# estim.int.effects() function
est.tmle <- estim.int.effects(interaction.ltmle, estimator = "tmle")
```

The results can be presented in a table following Know and VanderWeele
recommendations (2012). The out.table object contains the table and the
interaction.effects object contains the additive, multiplicative
interaction effects and the RERI.

The interaction table can be rendered using the `kableExtra` package.

``` r
table_inter <- out.int.table(int.r = est.tmle)

table_inter$out.table
table_inter$interaction.effects

library(kableExtra)
kbl(table_inter$out.table,
    caption = "Interaction effects estimated by TMLE") %>%
  kable_classic() %>%
  footnote(general = table_inter$interaction.effects)
```

We can also plot the results using the `out.int.fig()` function from the
output of the `estim.int.effects()` function.

``` r
out.int.fig(est.tmle)
```

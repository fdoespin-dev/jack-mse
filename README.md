# jack-mse <img src='https://www.ifop.cl/wp-content/themes/ifop-2024-07/assets/img/logo_ifop_transparente.png' align="right" height="44" />

<!-- badges: start -->

Development of operating models and empirical management procedures based on the acoustic survey of northern Chile for the jack mackerel fishery in the southeastern Pacific Ocean

The developed codes contain the main configurations designed to implement the base operating model with four fleets and seven abundance indices in [openMSE](https://github.com/Blue-Matter/openMSE) software, based on the current stock assessment model for jack mackerel, which assumes a single jack mackerel stock unit in the southeastern Pacific. In addition, three harvest control rules (HCRs) with two reference catch levels are coded for implementation in six HCR-based management procedures. These management procedures include a catch stabilization band that allows no more than a 20% increase or a 10% decrease in the catch.

## Example

Execute the following script to make the MOM in [R](https://www.r-project.org/),

```{r example, message=FALSE}
library(openMSE)
source("01-make-om-4-fleet.R")
```

This script creates the operating models from the stock assessment model, [JJM](https://github.com/SPRFMO/). Then, run the `03-simulate-mp.R`, this script calls the HCR-based management procedures `99-function-mmp.R` and run simulations.

```{r example, message=FALSE}
source("03-simulate-mp.R")
```

To generate the figures, run the following script:

```{r example, message=FALSE}
source("06-make-figures.R")
```

Projection plot,

<img src="/figuras/projection_SB.png" width="100%" />





# API

## Estimation

```@docs
panellp
PanelLocalProjection
WithinFE
SplitPanelJackknife
CDOLS
AlternatingProjections
OneSweep
```

## Formula helpers

```@docs
leads
ldiff
anchor
cumul
lags
firstdiff
fe
```

## Inference

```@docs
Robust
Clustered
DriscollKraaySE
vcov(::PanelLocalProjections.PanelVcov, ::PanelLocalProjection)
PanelLPVcov
stderror(::PanelLPVcov)
```

## Results

```@docs
summarize
IRFSummary
coefpath
nobspath
horizons
esample
irfplot
irfplot!
```

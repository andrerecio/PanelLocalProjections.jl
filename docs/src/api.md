# API

## Estimation

```@docs
panellp
PanelLocalProjection
WithinFE
SplitPanelJackknife
CDOLS
LagAugmented
ICLags
CappedLags
ConstantLags
LagSelection
lagorder
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
TLAHR
CosineHAR
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
synthetic
irfplot
irfplot!
```


selind <- tail(which(abscaug <= p),1)
xcap <- p - abscaug[selind]
if (f[selind] > -Inf) {
  exp(f[selind] + (f[selind + 1] - f[selind])/(abscaug[selind + 1] - abscaug[selind]) * xcap)
} else 0
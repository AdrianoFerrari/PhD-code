from autoscript import *

generate_jobs(
# name
  "stab1",
# conditions
  "(N == 20 and charge == 2.0) or (N == 40 and charge == 1.0)",

# Ns
  [20,40],

# Nls
  [20],

# charges
  [1.0,2.0],

# eps
  [1000.0],

# hs
  [1.0],

# hths
  [1.0],

# Lmaxs
  [0.4],

# T ----
  [21000],

# kTs
  [1.409],

# Rs
  [0.5],

# Po ----
  [100],

# Fo ----
  [10],

# seeds
  range(1,11),

# Dxs
  [0.01],

# Dxcs
  [0.01],

# turns
  [0.0],

# Amps
  [0.1],

# wavelengths
  [6.8/1,6.8/2,6.8/3,6.8/4,6.8/5,6.8/9,6.8/10,6.8/11,6.8/12],

# Lys
  [6.8],

# kcs
  [0],

# sigma particles
  [0.18],

# sigma chain-chain
  [0.18]
)

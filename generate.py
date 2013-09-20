from autoscript import *

generate_jobs(
# name
  "stab1",
# conditions
  "True",

# Ns
  [20],

# Nls
  [20],

# charges
  [2.0],

# eps
  [1000.0],

# hs
  [1.0],

# hths
  [1.0],

# Lmaxs
  [0.4],

# T ----
  [51000],

# kTs
  [1.409],

# Rs
  [0.5],

# Po ----
  [100],

# Fo ----
  [10],

# seeds
  range(1,6),

# Dxs
  [0.01],

# Dxcs
  [0.01],

# turns
  [0.0],

# Amps
  [0.0],

# wavelengths
  [1.0],

# Lys
  [6.8],

# kcs
  [0],

# sigma particles
  [0.18],

# sigma chain-chain
  [0.18],

# ebps
  [0.0,1.0,10.0,100.0,1000.0,10000.0]
)

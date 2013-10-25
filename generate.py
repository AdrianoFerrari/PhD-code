from autoscript import *

generate_jobs(
# name
  "helix6",
# conditions
  "True",

# Ns
  [20],

# Nls
  [10],

# charges
  [1.0],

# eps
  [1000.0],

# hs
  [1.0],

# hths
  [1.0],

# Lmaxs
  [0.4],

# T ----
  [161000],

# kTs
  [0.035],

# Rs
  [0.5],

# Po ----
  [1000],

# Fo ----
  [1000],

# seeds
  range(1,61),

# Dxs
  [0.008],

# Dxcs
  [0.01],

# turns
  [0.0],

# Amps
  [0.0],

# wavelengths
  [1.0],

# Lys
  [3.4],

# kcs
  [0],

# sigma particles
  [0.18],

# sigma chain-chain
  [0.18],

# ebps
  [0.0,0.01,0.025,0.05,0.1,0.25,0.5,1.0,2.5,5.0]
)

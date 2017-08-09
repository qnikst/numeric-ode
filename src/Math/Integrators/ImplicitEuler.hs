module Math.Integrators.ImplicitEuler
  ( implicitEuler
  ) where

import Linear

import Math.Integrators.Implicit

eps :: Floating a => a
eps = 1e-14

implicitEuler :: (Metric f, Ord a, Additive f, Num (f a), Floating a)
              => (f a -> f a) -> a -> f a -> f a
implicitEuler f = \h y ->
  fixedPoint (\x -> y ^+^ (h *^ (f x))) (\x1 x2 -> breakNormIR (x1^-^x2) eps) y

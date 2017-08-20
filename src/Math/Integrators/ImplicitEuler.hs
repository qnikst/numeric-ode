module Math.Integrators.ImplicitEuler
  ( implicitEuler
  ) where

import Linear

import Math.Integrators.Implicit

eps :: Floating a => a
eps = 1e-14

-- | Integrator of the form:
--
-- \[ y_{n+1} = y_n + h * f(\frac{y_n+y_{n+1}}{2}). \]
--
-- This is a symmetric method of order 1.
implicitEuler :: (Metric f, Ord a, Floating a)
              => (f a -> f a) -> a -> f a -> f a
implicitEuler f = \h y ->
  fixedPoint (\x -> y ^+^ (h *^ (f x))) (\x1 x2 -> breakNormIR (x1^-^x2) eps) y

module Math.Integrators.ImplicitEuler
  ( implicitEuler
  ) where

import Linear

import Math.Integrators.Implicit

-- | Integrator of the form:
--
-- \[ y_{n+1} = y_n + h * f(\frac{y_n+y_{n+1}}{2}). \]
--
-- This is a symmetric method of order 1.
implicitEuler :: (Metric f, Floating a)
              => (f a -> f a) -> Implicit f a
implicitEuler f = \h x0 x1 -> x0 ^+^ (h *^ (f x1))

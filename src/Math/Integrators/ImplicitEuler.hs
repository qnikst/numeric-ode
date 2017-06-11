{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.ImplicitEuler
  where

import Data.VectorSpace

import Math.Integrators.Implicit
import Math.Integrators.Internal

eps :: Double
eps = 1e-14::Double

implicitEuler :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> (a -> Double) -> Integrator a
implicitEuler f norm = \h y -> fixedPoint (\x -> y ^+^ (realToFrac h) *^ (f x)) (\x1 x2 -> breakNormR eps (norm (x1^-^x2))) y

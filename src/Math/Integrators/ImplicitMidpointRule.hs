{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.ImplicitMidpointRule
    where

import Data.Function
import Data.VectorSpace

import Math.Integrators.Implicit
import Math.Integrators.Internal

eps = 1e-14::Double

imr :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> (a -> Double) -> Integrator a
imr f norm = \h y -> fixedPoint (\x -> y ^+^ (realToFrac h) *^ ( f ( (y^+^x)^/2) )) (\x1 x2 -> breakNormR eps (norm (x1 ^-^ x2))) y

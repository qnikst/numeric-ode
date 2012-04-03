{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.ExplicitEuler
    where

import Math.Integrators.Internal
import Data.VectorSpace

explicitEuler :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> Integrator a
explicitEuler f = \h y -> y ^+^ ((realToFrac h) *^ (f y))

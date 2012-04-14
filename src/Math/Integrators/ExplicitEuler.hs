{-# LANGUAGE FlexibleContexts #-}
-- | Simplest integrations method used by Euler. 
-- properties: 
-- order:    1
--
module Math.Integrators.ExplicitEuler
    where

import Math.Integrators.Internal
import Data.VectorSpace

-- | Integrator of form
--  $$\Phi[h] : y_{n+1} = y_n + h f (y_n) $$
explicitEuler :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> Integrator a
explicitEuler f = \h y -> y ^+^ ((realToFrac h) *^ (f y))

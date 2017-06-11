{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module: Math.Integrators.ExplicitEuler.
--
-- Basic integrator using Euler method. It has following properies:
--   * allows to solve systems of the first order
--   * this method is not symplectic and tends to loose energy
--
module Math.Integrators.ExplicitEuler
    where

import Math.Integrators.Internal
import Data.VectorSpace

-- | Integrator of form
--
--  \[ \Phi[h] : y_{n+1} = y_n + h f (y_n) \]
explicitEuler :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> Integrator a
explicitEuler f = \h y -> y ^+^ ((realToFrac h) *^ (f y))

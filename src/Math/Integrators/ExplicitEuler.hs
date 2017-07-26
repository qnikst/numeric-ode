{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module: Math.Integrators.ExplicitEuler
--
-- Basic integrator using Euler method. It has following properies:
--
--   * Allows to solve systems of the first order.
--
--   * This method is not symplectic and tends to lose energy.
--
module Math.Integrators.ExplicitEuler
    (explicitEuler)
    where

import Linear

-- | Integrator of form: \( \Phi[h] : y_{n+1} = y_n + h f (y_n) \)
explicitEuler :: (Floating a, Additive f) => (f a -> f a) -> a -> f a -> f a
explicitEuler f = \h y -> y ^+^ h *^ (f y)

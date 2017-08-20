{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.ImplicitMidpointRule
  ( imr
  ) where

import Linear

import Math.Integrators.Implicit

imr :: (Metric f, Floating a)
    => (f a -> f a) -> Implicit f a
imr f = \h x0 x1 -> x0 ^+^ h *^ f((x0 ^+^ x1) ^/ 2)

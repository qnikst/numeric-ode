{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.ImplicitMidpointRule
  ( imr
  ) where

import Linear

import Math.Integrators.Implicit

eps :: Floating a => a
eps = 1e-14

imr :: (Metric f, Num (f a), Floating a, Ord a)
    => (f a -> f a) -> a -> f a -> f a
imr f = \h y ->
  fixedPoint (\x -> y ^+^ h *^ ( f ( (y^+^x)^/2) ))
             (\x1 x2 -> breakNormIR (x1 ^-^ x2) eps)
             y

{-# LANGUAGE FlexibleContexts #-}
-- | Helpers for implicit integration methods
--
-- This module is under heavy development and will likely to be
-- changed.
--
module Math.Integrators.Implicit
    ( -- * Types
      ImplicitSolver
      -- * Solvers
    , fixedPointSolver
    , fixedPoint
      -- * Helpers
    , breakNormR
    , breakNormIR
    )
    where

-- TODO Add possibility to make function to create initial value.
-- TODO add possibility to break on step
-- TODO: add possibility to add different initial value based
--          on y0, f
-- TODO: add seq-pseq to make this stuff strict
-- TODO: add Newton iterations

import Linear
import Control.Lens

-- | Implicit solver type
type ImplicitSolver a = (a -> a)                    -- ^ implicit method
                        -> (Int -> a -> a -> Bool)  -- ^ breakRule
                        -> a                        -- ^ initial value
                        -> a                        -- ^ final value

-- | Fixed point method it iterates function f until it will break "" will
-- be reached then it returns one but last iteration
--
fixedPointSolver :: ImplicitSolver a
fixedPointSolver f break' y0 = inner 0 y0
    where
        inner i y = let y' = f y
                        i' = i+1
                    in if break' i y y'
                           then y'
                           else inner i' y'

fixedPoint :: (a -> a)            -- ^ function
              -> (a -> a -> Bool) -- ^ break rule
              -> a                -- ^ initial value
              -> a                -- ^ result
fixedPoint f break' y0 =
    let y1 = f y0
    in if break' y0 y1
        then y0
        else fixedPoint f break' y1

-- | simple break rule that will break evaluatioin when value less then Eps
breakNormR :: Double -> Double -> Bool
breakNormR eps y =  abs y < eps

-- | same as @breakNormR@ but assume that inner type is an
-- instance of InnerField, so it's possible to use innerproduct to find norm
-- N.B function uses $||v||^2 < eps$, so epsilon should be pre evaluated
breakNormIR :: (Metric f, Floating a, Ord a, Num (f a)) => f a -> a -> Bool
breakNormIR v eps = quadrance v < eps




{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
-- | Helpers for implicit integration methods
--
-- Implicit methods are set of methods where next iteration
-- depends on the result of the next iteration.
--
-- Internals of the implicit method ('Implicit') takes an initial value and first
-- approximation to the value. User code should decide how to select first
-- approximation as complexity and efficiency of the method may depend on that
-- choice. Also user is free to select any strategy for running implicit method.
--
-- This module provides few basic methods for running iterations. But user have to
-- choose one that suits his needs best of all.
--
-- In order to run an implicit method one have to convert that to explicit one
-- that includes two steps:
--
--    1. Define how to find a required approximation
--    2. Define how to select initial step.
--
-- This module provides some helpers for that task.
module Math.Integrators.Implicit
    ( -- * Types
      Implicit
      -- * Solvers
    , fixed
    , nth
    , atNorm
      -- * Initial step
    , initExplicit
    )
    where

-- TODO Add possibility to make function to create initial value.
-- TODO add possibility to break on step
-- TODO: add possibility to add different initial value based
--          on y0, f
-- TODO: add seq-pseq to make this stuff strict
-- TODO: add Newton iterations

import Linear

-- | General helper type for the implicit method.
type Implicit f a
      = a -- ^ Delta. \(h\)
      -> f a  -- ^ Initial value \(x\)
      -> f a  -- ^ First approximation to the next value \(x'\)
      -> f a  -- ^ Result of application of the single step.


-- | Infinite list of fixed point iterations.
fixed :: Implicit f a -> a -> f a -> f a -> [f a] -- XXX: return stream instead, we are inifite?
fixed f dx x0 x1 = snd <$> iterate (\(x0',x1')  -> (x0, f dx x0' x1'))
                                   (x0,x1)

-- | Take \(n\)-th approximation.
nth :: Int -> Implicit f a -> a -> f a -> f a -> f a
nth n f dx x0 x1 = (fixed f dx x0 x1) !! n

atNorm :: forall f a . (Num (f a), Num a, Metric f, Ord a)
       => Implicit f a -> a -> a -> f a -> f a -> f a
atNorm f eps dx x0 x1 = undefined $ takeWhile check $ zip v (tail v)
  where
    v :: [f a]
    v = fixed f dx x0 x1
    check (a, b) = quadrance (a - b) >= eps

-- |
-- Helper that is using explicit method in order to find initial approximation.
-- In general case it's useful to use some basic method, for example euler
-- method.
initExplicit :: (a -> f a -> f a) -- ^ Explicit integrator.
         -> (a -> f a -> f a -> f a) -- ^ Implicit method solver.
         -> a                        -- ^ \(h\)
         -> f a                      -- ^ Initial value
         -> f a
initExplicit expl impl dx x0 = impl dx x0 (expl dx x0)

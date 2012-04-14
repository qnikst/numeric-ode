-- | Math integrators if a high level module for different ODE integrators
--   This module provides high-level wrappers over different integration methods
--
module Math.Integrators 
    where

import Data.Vector (Vector,(!))
import Data.Vector.Mutable
import Control.Monad.Primitive
import Control.Monad (liftM2)
import qualified Data.Vector as V

import Math.Integrators.Internal

{-|
 Integrate ODE equation using fixed steps set by a vector, and returns a vector
 of solutions corrensdonded to times that was requested.
 It takes Vector of time points as a parameter and returns a vector of results
 -}
integrateV :: PrimMonad m => Integrator a       -- ^ Internal integrator
                          -> a                  -- ^ initial  value
                          -> Vector Double      -- ^ vector of time points
                          -> m (Vector a)       -- ^ vector of solution
integrateV integrator initial times = do
    out <- new (V.length times) 
    write out 0 initial
    compute initial 1 out
    V.unsafeFreeze out
    where
        compute y i out | i == V.length times = return () 
                        | otherwise = do
            let h  = (times ! i) - (times ! (i-1))
                y' = integrator h y
            write out i y'
            compute y' (i+1) out


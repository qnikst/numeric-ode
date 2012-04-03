module Math.Integrators 
    where

import Data.Default
import Data.Vector (Vector,(!))
import Data.Vector.Mutable
import Control.Monad.Primitive
import Control.Monad (liftM2)
import qualified Data.Vector as V

import Math.Integrators.Internal
--type Integrator a = Double -> a -> a

integrateV :: PrimMonad m => Integrator a -> a -> Vector Double -> m (Vector a)
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

{-
integrateArray :: (a -> b) -> (Array a) -> a -> Integrator a -> Array (a,b)
integrateArray f ts x0 = do
-}   

{-# LANGUAGE RankNTypes #-}
module Math.Integrators.Internal
    where

--import Data.VectorSpace

-- | Integrator function
--
--    \Phi [h] |-> y_n -> y_{n+1}
-- 
type Integrator a = Double -> a -> a


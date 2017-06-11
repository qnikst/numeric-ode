{-# LANGUAGE FlexibleContexts #-}
-- | Stormer-Verlet has next properties
--
--      * the method is of order 2
--
--      * this is sympletic method
--
--      * method is symmetric
--
--      * ? for separate Hamiltonians T(p) + U(v) this method is explicit //not in code
--
--      * this method excatly conserves quadratic first integrals $p^T C$ (angular momentum of N-body system)
--


module Math.Integrators.StormerVerlet
    where

import Linear
import Control.Lens

-- | Stormer-Verlet integration scheme for system: 
--      
--      q'' = f(q)
--
stormerVerlet2 :: (Applicative f, Num (f a), Fractional a)
               => (f a -> f a)              -- ^ function f
               -> a
               -> V2 (f a)
               -> V2 (f a)
stormerVerlet2 f h prev =
    let h'  = h
        h2' = 0.5 *h
        p1  = prev ^. _x + pure h2' * (f (prev ^. _y))
        q'  = prev ^. _y + pure h' * p1
        p'  = p1 + pure h2' * (f q')
    in V2 p' q'


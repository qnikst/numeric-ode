{-# LANGUAGE FlexibleContexts #-}
-- Stormer-Verlet has next properties
-- * the method is of order 2
-- * this is sympletic method
-- * method is symmetric
-- * ? for separate Hamiltonians T(p) + U(v) this method is explicit //not in code
-- * this method excatly conserves quadratic first integrals $p^T C$ (angular momentum of N-body system)
--


module Math.Integrators.StormerVerlet
    where

import Data.AdditiveGroup
import Data.VectorSpace

import Math.Integrators.Internal
-- | St\"ormer-Verlet integration scheme for equatio q''=f(q)
--
-- p_{n+1/2} = p_{n-1/2} + h *f(q)
-- stormerVerlet1 :: (VectorSpace a, Floating (Scalar a)) -> (a->a) -> Integrator a
-- stormerVerlet1 = \h p -> -- undefined
    


-- | St\"ormer-Verlet integration scheme for system:
-- p'  = f(q)
-- q' = p
stormerVerlet2 :: (VectorSpace a, Floating (Scalar a)) => (a -> a) -> Integrator (a,a)
stormerVerlet2 f = \h (p,q) -> 
    let h'  = realToFrac h
        h2' = realToFrac $! 0.5 *h
        p1  = p ^+^ (h2') *^ (f q)
        q'  = q ^+^ h' *^ p1
        p'  = p1 ^+^ h2' *^ (f q')
    in (p',q')


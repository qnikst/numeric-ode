{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module: Math.Integrators.StormerVerlet
--
-- Stormer-Verlet has following properties
--
--  * the method is of order 2
--  * this is sympletic method
--  * method is symmetric
--  * this method excatly conserves quadratic first integrals \(p^T C\) (angular momentum of \(N\)-body system)
--
module Math.Integrators.StormerVerlet
    ( oneStepH98
    , stormerVerlet2
    ) where

import Linear
import Control.Lens

-- | Stormer-Verlet integration scheme for systems of the form
oneStepH98 :: (Applicative f, Num (f a), Fractional a) =>
              a            -- ^ Step size
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial q}\)
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial p}\)
           -> V2 (f a)     -- ^ Current \((p, q)\) as a 2-dimensional vector
           -> V2 (f a)     -- ^ New \((p, q)\) as a 2-dimensional vector
oneStepH98 hh nablaQ nablaP prev = V2 qNew pNew
  where
    h2   = hh / 2
    hhs  = pure hh
    hh2s = pure h2
    qsPrev = prev ^. _x
    psPrev = prev ^. _y
    pp2  = psPrev - hh2s * nablaQ qsPrev
    qNew = qsPrev + hhs * nablaP pp2
    pNew = pp2 - hh2s * nablaQ qNew

-- | Stormer-Verlet integration scheme for system: 
--      
--   \[   q'' = f(q) \]
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


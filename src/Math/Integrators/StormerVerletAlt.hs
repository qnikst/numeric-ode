{-# OPTIONS_GHC -Wall                   #-}
{-# OPTIONS_GHC -fno-warn-type-defaults #-}

module Math.Integrators.StormerVerletAlt where

import Linear
import Control.Lens

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

-- And now can apply this to the two body problem with the following
-- derivatives of the Hamiltonian.

nablaQ' :: V2 Double -> V2 Double
nablaQ' qs = V2 (qq1 / r) (qq2 / r)
  where
    qq1 = qs ^. _x
    qq2 = qs ^. _y
    r   = (qq1 ^ 2 + qq2 ^ 2) ** (3/2)

nablaP' :: V2 Double -> V2 Double
nablaP' ps = ps

e, q10, q20, p10, p20 :: Double
e = 0.6
q10 = 1 - e
q20 = 0.0
p10 = 0.0
p20 = sqrt ((1 + e) / (1 - e))

inits :: V2 (V2 Double)
inits = V2 (V2 q10 q20) (V2 p10 p20)

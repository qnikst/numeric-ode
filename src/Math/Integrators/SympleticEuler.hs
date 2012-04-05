{-# LANGUAGE FlexibleContexts #-}
module Math.Integrators.SympleticEuler
    where
import Data.Function
import Data.VectorSpace

import Math.Integrators.Implicit
import Math.Integrators.Internal
import Control.Parallel

eps :: Double
eps = 1e-10

sympleticEuler1 ::(VectorSpace a, Floating (Scalar a)) => ((a->a->a),(a->a->a)) -> (a -> Double) -> Integrator (a,a)
sympleticEuler1 (f,g) norm = \h (u,v) ->
    let u' = v' `pseq` u ^+^ ( (realToFrac h) *^ (f u v') )
        v' = fixedPoint (\x -> v ^+^ (realToFrac h) *^ (g u x)) (\x1 x2 -> breakNormR eps (norm (x1^-^x2))) v
    in (u',v')

{-
sEuler2 :: ((a->a->a),(a->a->a)) -> Double -> (a,a) -> (a,a)
sEuler2 (a,b) h (u,v) =
    let u' = u + ( h * (a u' v) )
        v' = v + ( h * (b u' v) )
    in (u',v')
-}

{-# LANGUAGE FlexibleContexts #-}
-- |
-- For partitioned systems
--
-- \[
--    \begin{array}{ccc}
--      \dot{u} & = & a(u,v) \\
--      \dot{v} & = & b(u,v)
--    \end{array}
-- \]
-- we can treat one variable by the implicit method and the other variable by the explicit
-- euler method.
module Math.Integrators.SympleticEuler
    where

import Linear
import Control.Lens

import Math.Integrators.Implicit

eps :: Floating a => a
eps = 1e-10


-- |
--
sympleticEuler1 :: (Metric f, Floating a, Ord a)
                => (f a -> f a -> f a)
                -> (f a -> f a -> f a)
                -> a                     -- ^ Step size
                -> V2 (f a)              -- ^ Current \((p,q)\) as a 2-dimentional vector
                -> V2 (f a)              -- ^ New \((p, q)\) as a 2-dimetional vector
sympleticEuler1 f g = \h prev ->
        -- explicit coordinate
    let u' = (prev^._x) ^+^  h *^ (f (prev^._x) v')
        -- implicit coordinate
        v' = fixedPoint (\x -> (prev^._y) ^+^ h *^ (g (prev^._x) x))
                        (\x1 x2 -> breakNormIR (x1^-^x2) eps) (prev^._x)
    in V2 u' v'

{-
sEuler2 :: ((a->a->a),(a->a->a)) -> Double -> (a,a) -> (a,a)
sEuler2 (a,b) h (u,v) =
    let u' = u + ( h * (a u' v) )
        v' = v + ( h * (b u' v) )
    in (u',v')
-}

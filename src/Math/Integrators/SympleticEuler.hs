{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
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
import Data.Functor.Compose
import Math.Integrators.Implicit
import Math.Integrators.ExplicitEuler
import Math.Integrators.ImplicitEuler

-- |
--
sympleticEuler1 :: forall f a  . (Metric f, Floating a)
                => (f a -> f a -> f a)
                -> (f a -> f a -> f a)
                -> Implicit (Compose V2 f) a
sympleticEuler1 a b = \h (Compose p) (Compose p') ->
  let u' = explicitEuler (\t -> a (p^._x) t) h (p^._y) :: f a
      v' = implicitEuler (\t -> b (p^._x) t) h (p^._y) (p'^._y) :: f a
  in Compose $ V2 u' v'

{-
sEuler2 :: ((a->a->a),(a->a->a)) -> Double -> (a,a) -> (a,a)
sEuler2 (a,b) h (u,v) =
    let u' = u + ( h * (a u' v) )
        v' = v + ( h * (b u' v) )
    in (u',v')
-}

{-# LANGUAGE FlexibleInstances, BangPatterns, TypeFamilies #-}
module Examples.LottkaVollterra
    where

import Control.Monad.ST
import Data.Vector (Vector,(!))
import qualified Data.Vector as V

import Math.Integrators
import Math.Integrators.ExplicitEuler

import Data.AdditiveGroup
import Data.VectorSpace

data T2 = T2 !Double !Double
    deriving (Eq,Show)
{-
   Lotka-Volterra model:

    \left\{ 
        \begin{array}{ccc}
            \dot{u} & = & u (v - 2) \\
            \dot{v} & = & v (1 - u) 
        \end{array}
     \right.
-}

-- | Equation in term of vectors
--lv :: Num a => (a,a) -> (a,a)

lv (T2 u v) = T2 (u * (v-2)) (v * (1-u))

-- | Initial conditions list
ics = [ V.fromList [2,2], V.fromList [4,8], V.fromList [4,2],V.fromList [6,2]]

-- | Initial conditions for testing
ic = T2 2 2

result = runST $ integrateV (explicitEuler lv) ic (V.enumFromStepN 0 0.12 100)

instance AdditiveGroup T2 where
    zeroV = T2 0 0
    negateV (T2 x y) = T2 (-x) (-y)
    (T2 a1 a2) ^+^ (T2 b1 b2) = T2 (a1+b1) (a2+b2)

instance VectorSpace T2 where
    type Scalar T2 = Double
    c *^ (T2 a1 a2) = T2 (c*a1) (c*a2)

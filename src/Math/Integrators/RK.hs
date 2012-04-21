{-# LANGUAGE QuasiQuotes, FlexibleContexts #-}
-- | Runge-Kutta module 
--   TODO: add description and history notes
--   add informations about methods properties
module Math.Integrators.RK
    ( -- * explicit methods
      rk45
    , rk46
      -- * implicit methods
    , gauss4
    , gauss6
    , lobattoIIIA4
    , lobattoIIIA6
    , lobattoIIIB4
    )
	where

import Math.Integrators.RK.Template
import Math.Integrators.RK.Types
import Math.Integrators.Internal
import Math.Integrators.Implicit
import Data.VectorSpace


rk45 :: (VectorSpace a, Floating (Scalar a)) => (Double -> a -> a) -> Integrator (Double,a)
rk45 = [qrk|
0   |
0.5 | 0.5
0.5 | 0   & 0.5
1   | 0   & 0   & 1
- - + - - - - - - - - 
    | 1/6 & 2/6 & 2/6 & 1/6
|]

rk46 :: (VectorSpace a, Floating (Scalar a)) => (Double -> a -> a) -> Integrator (Double,a)
rk46 = [qrk|
0   |
1/3 |  1/3
2/3 | -1/3 & 1
1   |  1   & -1  & 1
- - + - - - - - - - - 
    | 1/8  & 3/8 & 3/8 & 1/8
|]


gauss4 :: (VectorSpace a, Floating (Scalar a)) => (ImplicitRkType (a,a)) -> (Double -> a -> a) -> Integrator (Double,a)
gauss4 = [qrk|
0.5 - sqrt(3)/6 | 0.25 & 0.25 - sqrt(3)/6
0.5 + sqrt(3)/6 | 0.25 + sqrt(3)/6 & 1/4
- - - - - - - - + - - - - - - - - - - - 
                | 0.5     & 0.5
|]

gauss6 :: (VectorSpace a, Floating (Scalar a)) => (ImplicitRkType (a,a,a)) -> (Double -> a -> a) -> Integrator (Double,a)
gauss6 = [qrk|
0.5 - sqrt(15)/10 | 5/36 & 2/9 - sqrt(15)/15 & 5/36 - sqrt(15)/30
0.5               | 5/36 + sqrt(15)/24 & 2/9 & 5/36 - sqrt(15)/24
0.5 + sqrt(15)/10 | 5/36 + sqrt(15)/30 & 2/9+sqrt(15)/15 & 5/36
- - - - - - - - - + - - - - - - - - - - - - - - - - - - - - - - -
                  | 5/18 & 4/9 & 5/18
|]

lobattoIIIA4 :: (VectorSpace a, Floating (Scalar a)) => (ImplicitRkType (a,a,a)) -> (Double -> a -> a) -> Integrator (Double,a)
lobattoIIIA4 = [qrk|
0   |  0   &   0  & 0
0.5 | 5/24 & 1/3  & -1/24
1   | 1/6  & 2/3  & 1/6
- - + - - - - - - - - - - 
    | 1/6  & 2/3  & 1/6
|]

lobattoIIIA6 :: (VectorSpace a, Floating (Scalar a)) => (ImplicitRkType (a,a,a,a)) -> (Double -> a -> a) -> Integrator (Double,a)
lobattoIIIA6 = [qrk|
0               | 0 & 0 & 0 & 0
(5-sqrt(5))/10  | (11+sqrt(5))/120 & (25-sqrt(5))/120 & (25 - 13 *sqrt(5)/120) & (-1+sqrt(5))/120
(5+sqrt(5))/10  | (11-sqrt(5))/120 & (25+13*sqrt(5))/120 & (25+sqrt(5))/120    & (-1-sqrt(5))/120
1               | 1/12 & 5/12 &  5/12 & 1/12
- - - - - - - - + - - - -
                | 1/12 & 5/12 & 5/12 & 1/12
|]

lobattoIIIB4 :: (VectorSpace a, Floating (Scalar a)) => (ImplicitRkType (a,a,a)) -> (Double -> a -> a) -> Integrator (Double,a)
lobattoIIIB4 = [qrk|
0   | 1/6 & -1/6 & 0
0.5 | 1/6 & 1/3  & 0
1   | 1/6 & 5/6  & 0
- - + - - - - - - - - -
    | 1/6 & 2/3  & 1/6
|]

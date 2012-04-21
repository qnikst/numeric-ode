{-# LANGUAGE BangPatterns,TypeFamilies #-}
import Data.AdditiveGroup
import Data.VectorSpace
import qualified Data.Vector as V
import Control.Monad.ST

import Debug.Trace 
import Math.Integrators.StormerVerlet
import Math.Integrators.ExplicitEuler
import Math.Integrators.ImplicitEuler
import Math.Integrators.Internal
import Math.Integrators


import Data.Accessor
import Graphics.Rendering.Chart
import Data.Colour.Names
import Data.Colour



data T3 = T3 !Double !Double deriving (Eq,Show)
data T2 = T2 !T3 !T3 deriving (Eq,Show)

instance AdditiveGroup T3 where
    zeroV = T3 zeroV zeroV
    negateV (T3 x y) = T3 (negateV x) (negateV y)
    (T3 x1 y1) ^+^ (T3 x2 y2) = T3 (x1^+^x2) (y1^+^y2)

instance VectorSpace T3 where
    type Scalar T3 = Double
    c *^ (T3 x y) = T3 (c*^x) (c*^y)

instance InnerSpace T3 where
    (T3 x1 y1) <.> (T3 x2 y2) = x1*x2+y1*y2

instance AdditiveGroup T2 where
   zeroV = T2 zeroV zeroV
   negateV (T2 x y) = T2 (negateV x) (negateV y)
   (T2 a1 a2) ^+^ (T2 b1 b2) = T2 (a1^+^b1) (a2^+^b2)

instance VectorSpace T2 where
   type Scalar T2 = Double
   c *^ (T2 a1 a2) = T2 (c*^a1) (c*^a2)


kepler :: T2 -> T2
kepler (T2 q1 q2) = 
    let r  = q2 ^-^ q1          -- q2 - q1
        ri = r <.> r            -- ||q2-q1||^2
        rr = ri * (sqrt ri)
        q1' = r ^/ rr 
        q2' = negateV q1' 
    in T2 q1' q2'


k2p = \h (v,x) ->
    let v' = (implicitEuler (\_ -> kepler x) norm) h v
        x' = (implicitEuler (\_ -> v')  norm) h x
    in (v',x')


norm (T2 x y) = x<.>x + y<.>y

norm2 (a,b) = (norm a)*(norm a)+(norm b)*(norm b)


tm = V.enumFromStepN 0 0.001 100000
result1 = runST $ integrateV (stormerVerlet2 kepler) ((T2 (T3 1 0) (T3 0 1)),(T2 (T3 (-1) (-1)) (T3 1 1))) tm
result2 = runST $ integrateV (k2p) ((T2 (T3 1 0) (T3 0 1)),(T2 (T3 (-1) 0) (T3 1 0))) tm

plot1 r = renderableToPNGFile (mychart) 640 480 "all.png"
    where
        mychart = toRenderable layout 
        layout  = layout1_title ^= "bugs"
                $ layout1_plots ^= [Left (toPlot plot1),Left (toPlot plot2)]
                $ defaultLayout1
        plot1  = plot_lines_style .> line_color ^= opaque blue
               $ plot_lines_title  ^= "body-1"
               $ plot_lines_values ^= [map (\(_,T2 (T3 x y) _) -> (x,y))  (V.toList r)]
               $ defaultPlotLines
        plot2  = plot_lines_style .> line_color ^= opaque red
               $ plot_lines_title ^= "body-2"
               $ plot_lines_values ^= [map (\(_,T2 _ (T3 x y)) -> (x,y))  (V.toList r)]
               $ defaultPlotLines

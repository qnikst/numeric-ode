{-# OPTIONS_GHC -Wall         #-}

{-# LANGUAGE BangPatterns     #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}

module Main (main) where

import qualified Data.Vector as V
import Control.Monad.ST

import Math.Integrators.StormerVerlet
import Math.Integrators.ImplicitEuler
import Math.Integrators

import qualified Linear as L

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Cairo
import Data.Colour.Names
import Data.Colour
import Data.Default.Class

import Diagrams.Prelude ((.~))

kepler' :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
kepler' (L.V2 q1 q2) =
    let r  = q2 L.^-^ q1          -- q2 - q1
        ri = r `L.dot` r          -- ||q2-q1||^2
        rr = ri * (sqrt ri)
        q1' = r / pure rr
        q2' = negate q1'
    in L.V2 q1' q2'

k2p' :: L.V3 Double ->
        (L.V2 (L.V3 Double), L.V2 (L.V3 Double)) ->
        (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
k2p' = \h (v,x) ->
         let v' = (implicitEuler (\_ -> kepler' x)) h v
             x' = (implicitEuler (\_ -> v')) h x
         in (v',x')

tm :: V.Vector Double
tm = V.enumFromStepN 0 0.0001 1000000

u :: (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
u = ((L.V2 (L.V3 0.1 0 0) (L.V3 (-0.1) 0 0)), (L.V2 (L.V3 (-1) (-1) 0) (L.V3 1 1 0)))

u' :: L.V2 (L.V2 (L.V3 Double))
u' = L.V2 (L.V2 (L.V3 0.1 0 0) (L.V3 (-0.1) 0 0)) (L.V2 (L.V3 (-1) (-1) 0) (L.V3 1 1 0))

result1 :: V.Vector (L.V2 (L.V2 (L.V3 Double)))
result1 = runST $ integrateV (\h -> stormerVerlet2 kepler' (pure h)) u' tm

result2 :: V.Vector (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
result2 = runST $ integrateV (\h -> k2p' (pure h)) u tm

plotA :: V.Vector (L.V2 (L.V2 (L.V3 Double))) -> IO (PickFn ())
plotA r = renderableToFile def "../diagrams/integrators_example1_big.png" mychart
    where
        mychart :: Renderable ()
        mychart = toRenderable layout
        layout  = layout_title .~ "bugs"
                $ layout_plots .~ [toPlot plot1 ,toPlot plot2]
                $ def
        plot1  = plot_lines_style . line_color .~ opaque blue
               $ plot_lines_title  .~ "body-1"
               $ plot_lines_values .~ [map (\(L.V2 _ (L.V2 (L.V3 x y _z) _)) -> (x,y))  (V.toList r)]
               $ def
        plot2  = plot_lines_style . line_color .~ opaque red
               $ plot_lines_title .~ "body-2"
               $ plot_lines_values .~ [map (\(L.V2 _ (L.V2 _ (L.V3 x y _z))) -> (x,y))  (V.toList r)]
               $ def

main :: IO ()
main = do
  _ <- plotA result1
  putStrLn "Finished"

{-# OPTIONS_GHC -Wall         #-}

module Main (main) where

import qualified Data.Vector as V
import Control.Monad.ST

import Math.Integrators.StormerVerlet
import Math.Integrators.ImplicitEuler
import Math.Integrators

import qualified Linear as L

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Rasterific.CmdLine
import Codec.Picture.Gif

import Plots


myaxis :: Axis B V2 Double
myaxis = r2Axis &~ do
  linePlot' $ map unp2 $ take 200 morePts

myaxis' :: Int -> Axis B V2 Double
myaxis' n = r2Axis &~ do
  linePlot (map unp2 $ take n morePts) $ do
    plotColor .= black
  xMin .= Just (-6.00)
  xMax .= Just ( 0.00)
  yMin .= Just (-1.01)
  yMax .= Just ( 0.40)

kepler :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
kepler (L.V2 q1 q2) =
    let r  = q2 L.^-^ q1          -- q2 - q1
        ri = r `L.dot` r          -- ||q2-q1||^2
        rr = ri * (sqrt ri)
        q1' = r / pure rr
        q2' = negate q1'
    in L.V2 q1' q2'

k2p :: L.V3 Double ->
        (L.V2 (L.V3 Double), L.V2 (L.V3 Double)) ->
        (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
k2p = \h (v,x) ->
         let v' = (implicitEuler (\_ -> kepler x)) h v
             x' = (implicitEuler (\_ -> v')) h x
         in (v',x')

tm :: V.Vector Double
tm = V.enumFromStepN 0 0.0001 1000000

-- u :: (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
-- u = ((L.V2 (L.V3 1.2 0 0) (L.V3 (-0.9) 0 0)), (L.V2 (L.V3 (-1) (-1) 0) (L.V3 1 1 0)))

u' :: L.V2 (L.V2 (L.V3 Double))
u' = L.V2 (L.V2 (L.V3 0.1 0 0) (L.V3 (-0.5) 0 0)) (L.V2 (L.V3 (-1.5) (-1.0) 0) (L.V3 0.5 1.0 0))

result1 :: V.Vector (L.V2 (L.V2 (L.V3 Double)))
result1 = runST $ integrateV (\h -> stormerVerlet2 kepler (pure h)) u' tm

-- result2 :: V.Vector (L.V2 (L.V3 Double), L.V2 (L.V3 Double))
-- result2 = runST $ integrateV (\h -> k2p (pure h)) u tm

each' :: Int -> [a] -> [a]
each' n = map head . takeWhile (not . null) . iterate (drop n)

morePts :: [P2 Double]
morePts = map p2 $ each' 1000 $ map (\(L.V2 _ (L.V2 (L.V3 x y _z) _)) -> (x,y))  (V.toList result1)

displayHeader :: FilePath -> [(Diagram B, GifDelay)] -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , GifOpts {_dither = False, _noLooping = False, _loopRepeat = Just 10}
             )

displayHeader' :: FilePath -> Diagram B -> IO ()
displayHeader' fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  displayHeader' "other/plot-line.png" (renderAxis myaxis # bg white)
  displayHeader "other/anim-line.gif" $ zip (map (bg white . renderAxis . myaxis') [0..80]) (repeat 10)
  putStrLn "Finished"

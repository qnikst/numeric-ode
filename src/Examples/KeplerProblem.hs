{-# LANGUAGE NegativeLiterals #-}
{-# LANGUAGE TypeFamilies     #-}
{-# LANGUAGE FlexibleContexts #-}

{-# OPTIONS_GHC -Wall         #-}

module Main (main) where

import qualified Data.Vector as V
import Control.Monad.ST

import Math.Integrators.StormerVerlet
import Math.Integrators

import qualified Linear as L
import Linear.V
import Data.Maybe ( fromJust )

import Diagrams.Prelude
import Diagrams.Backend.CmdLine
import Diagrams.Backend.Rasterific.CmdLine

import Control.Monad
import Control.Monad.State.Class

import Plots


gConst :: Double
gConst = 6.67384e-11

nStepsTwoPlanets :: Int
nStepsTwoPlanets = 44

stepTwoPlanets :: Double
stepTwoPlanets = 24 * 60 * 60 * 100

sunMass, jupiterMass :: Double
sunMass     = 1.9889e30
jupiterMass = 1.8986e27

jupiterPerihelion :: Double
jupiterPerihelion = 7.405736e11

jupiterV :: [Double]
jupiterV = [-1.0965244901087316e02, -1.3710001990210707e04, 0.0]

jupiterQ :: [Double]
jupiterQ = [negate jupiterPerihelion, 0.0, 0.0]

sunV :: [Double]
sunV = [0.0, 0.0, 0.0]

sunQ :: [Double]
sunQ = [0.0, 0.0, 0.0]

tm :: V.Vector Double
tm = V.enumFromStepN 0 stepTwoPlanets nStepsTwoPlanets

kepler :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
kepler (L.V2 q1 q2) =
    let r  = q2 L.^-^ q1
        ri = r `L.dot` r
        rr = ri * (sqrt ri)
        q1' = pure gConst * r / pure rr
        q2' = negate q1'
        q1'' = q1' * pure sunMass
        q2'' = q2' * pure jupiterMass
    in L.V2 q1'' q2''

listToV3 :: [a] -> L.V3 a
listToV3 [x, y, z] = fromV . fromJust . fromVector . V.fromList $ [x, y, z]
listToV3 xs = error $ "Only supply 3 elements not: " ++ show (length xs)

initPQs :: L.V2 (L.V2 (L.V3 Double))
initPQs = L.V2 (L.V2 (listToV3 jupiterV) (listToV3 sunV))
                (L.V2 (listToV3 jupiterQ) (listToV3 sunQ))

result1 :: V.Vector (L.V2 (L.V2 (L.V3 Double)))
result1 = runST $ integrateV (\h -> stormerVerlet2 kepler (pure h)) initPQs tm

preMorePts :: [(Double, Double)]
preMorePts = map (\(L.V2 _ (L.V2 (L.V3 x y _z) _)) -> (x,y))  (V.toList result1)

morePts :: [P2 Double]
morePts = map p2 $ preMorePts

addPoint :: (Plotable (Diagram B) b, MonadState (Axis b V2 Double) m) =>
            Double -> (Double, Double) -> m ()
addPoint o (x, y) = addPlotable'
                    ((circle 1e11 :: Diagram B) #
                     fc brown #
                     opacity o #
                     translate (r2 (x, y)))

jSaxis :: Axis B V2 Double
jSaxis = r2Axis &~ do
  addPlotable' ((circle 1e11 :: Diagram B) # fc yellow)
  let l = length preMorePts
  let os = [0.05,0.1..]
  let ps = take (l `div` 4) [0,4..]
  zipWithM_ addPoint os (map (preMorePts!!) ps)
  linePlot' $ map unp2 $ take 200 morePts

displayHeader :: FilePath -> Diagram B -> IO ()
displayHeader fn =
  mainRender ( DiagramOpts (Just 900) (Just 700) fn
             , DiagramLoopOpts False Nothing 0
             )

main :: IO ()
main = do
  displayHeader "other/jupiter-sun-line.png" (renderAxis jSaxis # bg ivory)
  putStrLn "Finished"


{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module: Math.Integrators.StormerVerlet
--
-- Störmer-Verlet is
--
--  * Of order 2
--  * A sympletic method
--  * Symmetric
--  * Excatly conserves quadratic first integrals \(p^T C\) (angular momentum of \(N\)-body system)
--
module Math.Integrators.StormerVerlet
    ( oneStepH98
    , stormerVerlet2
    ) where

import Linear
import Control.Lens

-- | Stormer-Verlet integration scheme for systems of the form
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

-- | Stormer-Verlet integration scheme for system: \(\ddot{\mathbf{q}} = f(\mathbf{q})\)
--
-- > {-# LANGUAGE NegativeLiterals      #-}
-- > {-# LANGUAGE TypeFamilies          #-}
-- > {-# LANGUAGE FlexibleContexts      #-}
-- > {-# LANGUAGE MultiParamTypeClasses #-}
-- > 
-- > import qualified Data.Vector as V
-- > import Control.Monad.ST
-- > 
-- > import Math.Integrators.StormerVerlet
-- > import Math.Integrators
-- > 
-- > import qualified Linear as L
-- > import Linear.V
-- > import Data.Maybe ( fromJust )
-- > 
-- > import Diagrams.Prelude
-- > 
-- > import Control.Monad
-- > import Control.Monad.State.Class
-- > 
-- > import Plots
-- > 
-- > 
-- > gConst :: Double
-- > gConst = 6.67384e-11
-- > 
-- > nStepsTwoPlanets :: Int
-- > nStepsTwoPlanets = 44
-- > 
-- > stepTwoPlanets :: Double
-- > stepTwoPlanets = 24 * 60 * 60 * 100
-- > 
-- > sunMass, jupiterMass :: Double
-- > sunMass     = 1.9889e30
-- > jupiterMass = 1.8986e27
-- > 
-- > jupiterPerihelion :: Double
-- > jupiterPerihelion = 7.405736e11
-- > 
-- > jupiterV :: [Double]
-- > jupiterV = [-1.0965244901087316e02, -1.3710001990210707e04, 0.0]
-- > 
-- > jupiterQ :: [Double]
-- > jupiterQ = [negate jupiterPerihelion, 0.0, 0.0]
-- > 
-- > sunV :: [Double]
-- > sunV = [0.0, 0.0, 0.0]
-- > 
-- > sunQ :: [Double]
-- > sunQ = [0.0, 0.0, 0.0]
-- > 
-- > tm :: V.Vector Double
-- > tm = V.enumFromStepN 0 stepTwoPlanets nStepsTwoPlanets
-- > 
-- > kepler :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
-- > kepler (L.V2 q1 q2) =
-- >     let r  = q2 L.^-^ q1
-- >         ri = r `L.dot` r
-- >         rr = ri * (sqrt ri)
-- >         q1' = pure gConst * r / pure rr
-- >         q2' = negate q1'
-- >         q1'' = q1' * pure sunMass
-- >         q2'' = q2' * pure jupiterMass
-- >     in L.V2 q1'' q2''
-- > 
-- > listToV3 :: [a] -> L.V3 a
-- > listToV3 [x, y, z] = fromV . fromJust . fromVector . V.fromList $ [x, y, z]
-- > listToV3 xs = error $ "Only supply 3 elements not: " ++ show (length xs)
-- > 
-- > initPQs :: L.V2 (L.V2 (L.V3 Double))
-- > initPQs = L.V2 (L.V2 (listToV3 jupiterV) (listToV3 sunV))
-- >                 (L.V2 (listToV3 jupiterQ) (listToV3 sunQ))
-- > 
-- > result1 :: V.Vector (L.V2 (L.V2 (L.V3 Double)))
-- > result1 = runST $ integrateV (\h -> stormerVerlet2 kepler (pure h)) initPQs tm
-- > 
-- > preMorePts :: [(Double, Double)]
-- > preMorePts = map (\(L.V2 _ (L.V2 (L.V3 x y _z) _)) -> (x,y))  (V.toList result1)
-- > 
-- > morePts :: [P2 Double]
-- > morePts = map p2 $ preMorePts
-- > 
-- > addPoint :: (Plotable (Diagram B) b, MonadState (Axis b V2 Double) m) =>
-- >             Double -> (Double, Double) -> m ()
-- > addPoint o (x, y) = addPlotable'
-- >                     ((circle 1e11 :: Diagram B) #
-- >                      fc brown #
-- >                      opacity o #
-- >                      translate (r2 (x, y)))
-- > 
-- > jSaxis :: Axis B V2 Double
-- > jSaxis = r2Axis &~ do
-- >   addPlotable' ((circle 1e11 :: Diagram B) # fc yellow)
-- >   let l = length preMorePts
-- >   let os = [0.05,0.1..]
-- >   let ps = take (l `div` 4) [0,4..]
-- >   zipWithM_ addPoint os (map (preMorePts!!) ps)
-- >   linePlot' $ map unp2 $ take 200 morePts
-- > 
-- > jupiterOrbit = renderAxis jSaxis # bg ivory
--
-- <<diagrams/src_Math_Integrators_StormerVerlet_jupiterOrbit.svg#diagram=jupiterOrbit&height=300&width=200>>
--
-- For example, for the \(n\)-body problem, the Hamiltonian is
--
-- \[
-- {\mathbb H} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} - \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
-- \]
--
-- Apply Hamilton's equations will give us \(2n\) first order
-- equations. To use 'stormerVerlet2' we need \(n\) second order
-- equations. In this case, the Lagrangian is easy
--
-- \[
-- {\mathcal{L}} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} + \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
-- \]
--
-- Applying Lagrange's equation
--
-- \[
-- \frac{\mathrm{d}}{\mathrm{d}t}\bigg(\frac{\partial{\mathcal{L}}}{\partial\dot{q}_j}\bigg) = \frac{\partial{\mathcal{L}}}{\partial{q}_j}
-- \]
--
-- gives
--
-- \[
-- m_j\ddot{q}_j = G\sum_{k \neq j}m_k m_j \frac{q_k - q_j}{\|q_k - q_j\|^3}
-- \]
--
-- For \(n = 2\) this gives
--
-- \[
-- \begin{aligned}
-- m_1\ddot{q}_1 &= G\frac{q_1 - q_2}{\|q_1 - q_2\|^3} \\
-- m_2\ddot{q}_2 &= G\frac{q_2 - q_1}{\|q_2 - q_1\|^3}
-- \end{aligned}
-- \]
--
-- > kepler' :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
-- > kepler' (L.V2 q1 q2) =
-- >     let r  = q2 L.^-^ q1          -- q2 - q1
-- >         ri = r `L.dot` r          -- ||q2-q1||^2
-- >         rr = ri * (sqrt ri)
-- >         q1' = r / pure rr
-- >         q2' = negate q1'
-- >     in L.V2 q1' q2'
--
stormerVerlet2 :: (Applicative f, Num (f a), Fractional a)
               => (f a -> f a)              -- ^ function f
               -> a
               -> V2 (f a)
               -> V2 (f a)
stormerVerlet2 f h prev =
    let h'  = h
        h2' = 0.5 *h
        p1  = prev ^. _x + pure h2' * (f (prev ^. _y))
        q'  = prev ^. _y + pure h' * p1
        p'  = p1 + pure h2' * (f q')
    in V2 p' q'


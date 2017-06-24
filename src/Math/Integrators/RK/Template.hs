{-# LANGUAGE TemplateHaskell #-}
{-# OPTIONS_GHC -Wwarn #-} -- this module will be removed in future versions


-- Let \((b_i)\),\((a_{ij},(i,j=1,\ldots,s)\)\) be real numbers and let
-- \((c=\sum_{j=1}a_{ij})\). An $s$-stage Runge-Kutta method is given by
--
--   \(\(k_i = f(t_0+ c_ih, y_0+h\sum_{j=1}^s a_{ij} k_j), i=1,\ldots,s\)\)
--   \(\(y_1 = y_0 + h\sum_{i=1}^s b_ik_i\)\)
--
-- Te coefficients are usually displayed as follows:
--
--    c_1 | a_11 .. a_{1s}
--    ... | ...     ...
--    c_s | a_s1 ... a_{ss}
--    ----+----------------
--        | b_1 ...  b_s
--
-- This module provide a template haskell based approach for method generation
-- all elements takes a table as a list:
--
--    [[c1, a_11 ... a_1s]
--    , ...
--    ,[c_s, a_s1 ... a_ss]
--    ,    [ b_1 ...  b_s]
--    ]
--
module Math.Integrators.RK.Template
    where

import Language.Haskell.TH.Syntax
import Language.Haskell.TH

newVar :: String -> Q (Name, ExpQ, PatQ)
newVar s = do
  z <- newName s
  return (z,varE z, varP z)

-- | Generate Runge-Kutta method based on the table.
generateRK :: Lift a => [[a]] -> ExpQ
generateRK matrix = do
    (_,fE,fP) <- newVar "f"
    (_,hE,hP) <- newVar "h"
    (_,y0E, y0P) <- newVar "y0"
    let makeKas = zipWith go [(0::Int)..] ks
          where
            go i name = funD name [clause [] (normalB (k i)) []]
            k i = [| $fE (t0 + $(mkC (cs !! i)) * $hE) ($y0E + $hE * $(sumOp (as !! i))) |]
        y1 = [| $y0E + $hE * $(sumOp bs) |]
    [| \ $fP $hP $y0P -> $(letE makeKas  [| $y1 |] )
     |]
  where
    -- split the matrix into parameters
    as = map (drop 1) $ init matrix
    bs = last matrix
    cs = concat $ take 1 $ init matrix
    s  = length matrix - 1
    ks = map (\i -> mkName ("k" ++ show i)) [0..(s-1)]
    -- Generate $k_i$ according to formula above
    mkC c = [| c |]
    sumOp aa = internal 0 aa
      where internal i [x]    = [| x * $(varE $ ks !! i)|]
            internal i (x:xs) = [| x * $(varE $ ks !! i) + $(internal (i+1) xs) |]
            internal _ []     = [| 0 |]


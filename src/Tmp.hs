module Tmp
	where

import Control.Parallel

data SumE = SumE Double Double
	deriving (Show, Eq, Read)

add :: SumE -> Double -> SumE
add (SumE y e) d = let e'= e + d
		       y'= y+e'
		       e''= e'+(y-y')
		   in y' `pseq` e'' `pseq` (SumE y' e'')

pty :: Double -> SumE
pty x = SumE x 0

out :: SumE -> Double
out (SumE x _) = x

module Math.Integrators.Implicit
    where


breakNormR :: Double -> Double -> Bool
breakNormR eps y =  y < eps

fixedPoint :: (a -> a) -> (a -> a -> Bool) -> a -> a
fixedPoint f break y0 = 
    let y1 = f y0
    in if break y0 y1
        then y0
        else fixedPoint f break y1

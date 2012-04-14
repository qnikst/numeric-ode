-- | Helpers for implicit integration methods
--
-- TODO: add possibility to break on step
-- TODO: add possibility to add different initial value based
--          on y0, f
-- TODO: add seq-pseq to make this stuff strict
-- TODO: add Newton iterations
module Math.Integrators.Implicit
    where

-- | simple break rule that will break evaluatioin when value less then Eps
breakNormR :: Double -> Double -> Bool
breakNormR eps y =  y < eps

-- | same as @breakNormR@ but assume that inner type is an 
-- instance of InnerField, so it's possible to use innerproduct to find norm
-- N.B function uses $||v||^2 < eps$, so epsilon should be pre evaluated
breakNormIR :: (InnerSpace a, VectorSpace (Floating b)) => a -> b -> Bool
breakNormIR v eps = (v <.> v) < eps


-- | Fixed point method it iterates function f until it will break "" will
-- be reached then it returns one but last iteration
--
fixedPoint :: (a -> a)            -- ^ function
              -> (a -> a -> Bool) -- ^ break rule
              -> a                -- ^ initial value
              -> a                -- ^ result
fixedPoint f break y0 = 
    let y1 = f y0
    in if break y0 y1
        then y0
        else fixedPoint f break y1

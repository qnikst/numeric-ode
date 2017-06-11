module Math.Integrators.RK.Internal
    ( MExp(..)
    , isExplicit
    )
    where


-- | Internal type that users by
data MExp = Delimeter | Row (Maybe Double,[Double]) deriving (Show,Eq)

isExplicit :: [MExp] -> Bool
isExplicit  = (all (\(i,(Row (_,x))) -> i> length x)) . (zip [1..]) . top
    where 
        top = takeWhile (/= Delimeter) 

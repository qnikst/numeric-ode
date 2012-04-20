{-# LANGUAGE TemplateHaskell #-}
module Math.Integrators.RK.Template
    where

import Data.Maybe

import Language.Haskell.TH.Quote
import Language.Haskell.TH

import Control.Monad

import Math.Integrators.RK.Types
import Math.Integrators.RK.Internal
import Math.Integrators.RK.Parser


qrk  :: QuasiQuoter
qrk  =  QuasiQuoter {quoteExp = x}
    where
        x s = rk $! readMatrixTable s

-----------------------------------------------------------
-- List of helpers
jv    = Just . VarE
jv'   = Just . varE
jld   = Just . LitE. {-DoublePrimL-} RationalL . toRational
ld    = LitE . RationalL . toRational
plus  = VarE $! mkName "+"
vplus = VarE $! mkName "^+^"
vmult = VarE $! mkName "*^"
mult  = VarE $! mkName "*"
ld'   = litE . RationalL . toRational
plus' = varE $! mkName "+"
vplus'= varE $! mkName "^+^"
vmult'= varE $! mkName "*^"
mult' = varE $! mkName "*"
vN s = varE (mkName s)

foldOp op = foldl1 (\x y -> infixE (Just x) op (Just y))

realToFracN = varE (mkName "realToFrac")
f = mkName "f"
t = mkName "t"
h = mkName "h"
y = mkName "y"

rk :: [MExp] -> Q Exp
rk mExp = do
    let (ab,_:c:[]) = break (==Delimeter) mExp
        lenA     = length ab
    let f = mkName "f"
        h = mkName "h"
        t = mkName "t"
        y = mkName "y"
    kn <- forM [1..lenA] (\_ -> newName "k")
    let kvv = zip kn ab
    ks' <- forM kvv $ \(k,r) -> do
            t <- rkt1 r kn
            return $ ValD (VarP k) (NormalB t) []
    y' <- rkt2 c kn
    return $ LamE [VarP f, VarP h, TupP [VarP t,VarP y]] $ LetE ks' y'
    where 
        rkt1 (Row (Just c,ls)) ks = do
            let f = mkName "f"
                h = mkName "h"
                t = mkName "t"
                y = mkName "y"
            let ft = InfixE (jv t) plus (Just $ InfixE (jld c) mult (jv h))
                st = if null ls 
                        then (VarE y)
                        else
                            InfixE (jv y) vplus
                                (Just $ InfixE (Just $ AppE (VarE $ mkName "realToFrac") (VarE h)) (vmult)
                                    (Just $ foldl1 (\x y -> InfixE (Just x) vplus (Just y)) $
                                                zipWith (\k l -> InfixE (Just $ AppE (VarE $ mkName "realToFrac") (ld l)) vmult (jv k)) ks ls
                                    )
                                )
            return $ AppE ( AppE (VarE f) ft) st
        rkt2 (Row (_,ls)) ks = do
            let f = mkName "f"
                y = mkName "y"
                h = mkName "h"
                t = mkName "t"
            return $ TupE 
                        [ InfixE (jv t) plus (jv h)
                        , InfixE (jv y) vplus
                          (Just $ InfixE (Just $ AppE (VarE $ mkName "realToFrac") (VarE h)) vmult
                            (Just $ foldl1 (\x y -> InfixE (Just x) vplus (Just y)) $
                                        zipWith (\k l -> InfixE (Just $ AppE (VarE $ mkName "realToFrac") (ld l)) vmult (jv k)) ks ls
                            )
                          )
                        ]


test = [Row (Just 1,[2,3]),Row (Just 4,[5,6]),Delimeter, Row (Nothing, [7,8])]

irk mExp = do
    let lena     = length "ab"
        tpy      = mkName "tpy"
        f        = mkName "f"
        h        = mkName "h"
        t        = mkName "t"
        y        = mkName "y"
    fpoint' <- fpoint mExp
    lamE [varP tpy, varP f, varP h, tupP [varP t,varP y]] $ 
        caseE (varE tpy) [match (varP $ mkName "Math.Integrators.RK.Types.FixedPoint breakRule") 
                                (normalB (fpointRun mExp))  
                                fpoint'
                         ,match (varP $ mkName "Math.Integrators.RK.Types.NewtonIteration") (normalB (varE $ mkName "undefined")) []
                         ]

fpointRun mExp = do
    let (ab,_:(Row (_,ls)):[]) = break (==Delimeter) mExp
        lenA     = length ab
    zs <- forM [1..lenA] (\_ -> newName "z")
    letE [valD   (tupP $ map varP zs)
                 (normalB $ 
                     appE 
                        (appE 
                            (appE 
                                (varE $! mkName "Math.Integrators.Implicit.fixedPointSolver") 
                                (varE $! mkName "method")
                            ) 
                            (varE $ mkName "breakRule")
                        ) 
                        (tupE $ map (litE . RationalL) $ replicate lenA 0) {- TODO: give avaliability to user -}
                     ) 
                  []]
         (appE (varE (mkName "solution")) (tupE $ map varE zs))

fpoint mExp = do
    let (ab,_:(Row (_,ls)):[]) = break (==Delimeter) mExp
        lenA     = length ab
    zs <- forM [1..lenA] (\_ -> newName "z")

    return $ 
        [ funD (mkName "method") [clause [varP y] (normalB $ letE (map (topRow zs) $! zip zs ab) (tupE $ map varE zs)) [] ]
        , funD (mkName "solution") [clause [tupP $ map varP zs] (normalB $ solutionRow zs ls) []]
        ]
    where
        topRow zs (x,(Row (Just c,ls))) = 
            valD (varP x) 
                 (normalB $ infixE (jv' h) 
                         vmult'
                         (Just $ foldOp vplus' $
                                    zipWith (\z l -> infixE (Just $ appE realToFracN (ld' l))
                                                            vmult'
                                                            (Just $ appE (appE (varE f) 
                                                                             (infixE (jv' t) 
                                                                                     plus' 
                                                                                     (Just $ infixE (jv' h) mult' (Just $ ld' c))
                                                                              )
                                                                          )
                                                                          (infixE (jv' y) vplus' (jv' z))
                                                            )
                                            )
                                            zs
                                            ls
                        )
                 )
                []
        solutionRow zs ls = 
                (tupE [infixE (jv' t) plus' (jv' h)
                      ,infixE (jv' y) 
                            vplus' 
                            (Just $ infixE (jv' h) 
                                vmult'
                                (Just $ foldOp vplus' 
                                          (zipWith (\z b -> infixE (Just $ appE realToFracN 
                                                                                (ld' b))
                                                                   vmult'
                                                                   (jv' z))
                                                   zs
                                                   ls
                                          ) 
                                )
                            )
                      ]
                  )
        {-
        
        method = do 
            let (ab,_:(Row (_,ls)):[]) = break (==Delimeter) mExp
                lenA     = length ab
            zs <- forM [1..lenA] (\_ -> newName "z")
            letE 
                [ funD (mkName "mymethod") 
                , valD (tupP $ map varP zs) 
                       (normalB $ 
                        appE 
                            (appE 
                                (appE 
                                    (varE $! mkName "Math.Integrators.Implicit.fixedPointSolver") 
                                    (appE (varE $! mkName "mymethod") (varE t))
                                ) 
                                (varE $ mkName "breakRule")
                            ) 
                            (tupE $ map (litE . RationalL) $ replicate lenA 0) {- TODO: give avaliability to user -}
                        ) 
                        [] 
                ] (tupE [infixE (jv' t) plus' (jv' h)
                        ,infixE (jv' y) 
                                vplus' 
                                (Just $ infixE (jv' h) 
                                        vmult'
                                        (Just $ foldl1 (\x y -> infixE (Just x) vplus' (Just y)) 
                                                       (zipWith (\z b -> infixE (Just $ appE (varE $ mkName "realToFrac")
                                                                                             (ld' b))
                                                                                vmult'
                                                                                (jv' z))
                                                                zs
                                                                ls
                                                        )
                                        )
                                )
                        ]
                  )
        -}
    {-
     - Z_in = h * a_j f^j (t_j, k_j)
     -}
{-
implicit = \t f implicit (t,y) -> 
    let (z1,z2,...zN) = implicit
                            (\(z1,z2..zN) ->
                                let z1' = if (z1 > eps) then h * sum (a_{ij} f (y+z1)) else z1
                                    z2' = if (z2 > eps) then h * sum (................) else z2
                                    ....
                            in tuple (z1',z2',z3'..zN'))
                            (z1,z2,z3,...,zN')
    in y ^+^ h *^ sum ( b_1 * (y ^+^ z1))
-}


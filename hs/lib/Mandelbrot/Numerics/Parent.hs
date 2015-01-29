{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Parent
  ( Parent(..)
  , parent
  , parents
  ) where

import Data.Ratio
  ( (%)
  )
import Data.Strict.Tuple
  ( Pair(..)
  )

import Mandelbrot.Numerics.Complex
import Mandelbrot.Numerics.Interior
import Mandelbrot.Numerics.Nucleus
import Mandelbrot.Numerics.Shape
import Mandelbrot.Numerics.Wucleus

data Parent r
  = NoParent !(Complex r)
  | Parent !(Complex r) !Rational !Int !(Complex r)
  deriving (Eq, Read, Show)

parent
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> Maybe (Parent r)
{-# SPECIALIZE parent :: Int -> Complex Double -> Maybe (Parent Double) #-}
parent p c0 = case shape p c0 of
  Nothing -> Nothing
  Just Cardioid -> case converge 2 64 $ interior p 1 c0 c0 of
    Nothing -> Nothing
    Just (_:!:r) -> Just (NoParent r)
  Just Circle -> case converge 2 64 $ interior p ((k-3)/k) c0 c0 of
    Just (z1:!:c1) -> case converge 2 64 $ interior p ((k-1)/k) z1 c1 of
      Just (_:!:c2) -> go (double c2 - c1) 1 (1/0) 0
      _ -> Nothing
    _ -> Nothing
    where
      k = 32
      go c1 q zq z
        | q >= p = Nothing
        | zq' < zq && p `mod` q == 0 = case converge 2 64 $ wucleus q c1 z' of
            Just z1 -> case zdz q c1 z1 1 of
              (_:!:dz1) | magnitudeSquared dz1 <= 1 ->
                let den = p `div` q
                    num = round (fromIntegral den * phase dz1 / (2 * pi)) `mod` den
                    t = toInteger num % toInteger den
                    i = cis (2 * pi * fromRational t)
                in  case converge 2 64 $ nucleus q c1 of
                      Just c2 -> case converge 2 64 $ interior q i c2 c2 of
                        Just (_:!:c3) -> Just $ Parent c3 t q c2
                        _ -> Nothing
                      _ -> Nothing
              _ -> go c1 (q + 1) zq' z'
            _ -> go c1 (q + 1) zq' z'
        | zq' < zq = go c1 (q + 1) zq' z'
        | otherwise = go c1 (q + 1) zq z'
        where
          z' = sqr z + c1
          zq' = magnitudeSquared z'
  where
    zdz q c !z !dz
      | q == 0 = z :!: dz
      | otherwise = zdz (q - 1) c (sqr z + c) (double (z * dz))

parents
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> [Parent r]
{-# SPECIALIZE parents :: Int -> Complex Double -> [Parent Double] #-}
parents p c = case parent p c of
  Nothing -> []
  Just q@(NoParent _) -> [q]
  Just q@(Parent _ _ p' c') -> q : parents p' c'

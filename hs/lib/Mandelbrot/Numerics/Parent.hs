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
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Maybe (Parent r)
parent p c0 = case shape p c0 of
  Nothing -> Nothing
  Just Cardioid -> case interior' p 1 c0 c0 of
    Nothing -> Nothing
    Just (_:!:r) -> Just (NoParent r)
  Just Circle -> case internalRay [] (1/128) (1/64) c0 c0 of
    Just ((_:!:root1):(_:!:root0):_) -> go (double root1 - root0) 1 (1/0) 0
    _ -> Nothing
    where
      go c1 q zq z
        | q >= p = Nothing
        | zq' < zq && p `mod` q == 0 = case wucleus' q c1 z' of
            Just z1 -> case zdz q c1 z1 1 of
              (_:!:dz1) | magnitudeSquared dz1 <= 1 ->
                let den = p `div` q
                    num = round (fromIntegral den * phase dz1 / (2 * pi)) `mod` den
                    t = toInteger num % toInteger den
                    i = cis (2 * pi * fromRational t)
                in  case nucleus' q c1 of
                      Just c2 -> case interior' q i c2 c2 of
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
    internalRay acc r inc !z !c
      | magnitudeSquared r >= 1 = Just acc
      | otherwise = case interior' p r z c of
          Just zc@(z':!:c') -> internalRay (zc:acc) (r + inc) inc z' c'
          _ -> Nothing
    nucleus' q c = converge $ nucleus q c
    wucleus' q c z = converge $ wucleus q c z
    interior' q i z c = converge2 $ interior q i z c
    converge = lastMay . filter finiteC . take 64
    converge2 = lastMay . filter finiteC2 . take 64
    finiteC2 (a:!:b) = finiteC a && finiteC b
    lastMay [] = Nothing
    lastMay xs = Just (last xs)
    zdz q c !z !dz
      | q == 0 = z :!: dz
      | otherwise = zdz (q - 1) c (sqr z + c) (double (z * dz))

parents
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> [Parent r]
parents p c = case parent p c of
  Nothing -> []
  Just q@(NoParent _) -> [q]
  Just q@(Parent _ _ p' c') -> q : parents p' c'

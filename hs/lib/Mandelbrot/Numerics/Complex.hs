{- |
Complex number operations useful in Mandelbrot set numeric algorithms.
-}
module Mandelbrot.Numerics.Complex
  ( module Data.Complex
  , Approx(..)
  , converge
  , Square(..)
  , magnitudeSquared
  , loss
  ) where

import Data.Complex
import Data.Strict.Tuple
  ( Pair(..)
  )
import Safe
  ( lastMay
  )

class Approx t where
  finite :: t -> Bool
  finite = not . notFinite
  notFinite :: t -> Bool
  notFinite = not . finite
  bashZero :: t -> t
  approxEq :: Int -> t -> t -> Bool

instance Approx Double where
  notFinite x = isNaN x || isInfinite x
  bashZero x
    | 1 + abs x == 1 = 0
    | otherwise = x
  approxEq b u v
    | x == 0 && y == 0 = True
    | x >  0 && y >  0 = loss   x    y  > k
    | x <  0 && y <  0 = loss (-x) (-y) > k
    | otherwise = False
    where
      x = bashZero u
      y = bashZero v
      k = fromIntegral (floatDigits x - b)

loss :: Double -> Double -> Double
loss p q = negate . logBase 2 $ 1 - min p q / max p q

instance Approx t => Approx (Complex t) where
  {-# SPECIALIZE instance Approx (Complex Double) #-}
  finite (x :+ y) = finite x && finite y
  bashZero (x :+ y) = bashZero x :+ bashZero y
  approxEq b (ux :+ uy) (vx :+ vy) = approxEq b ux uy && approxEq b vx vy

instance (Approx s, Approx t) => Approx (Pair s t) where
  {-# SPECIALIZE instance Approx (Pair (Complex Double) (Complex Double)) #-}
  finite (a :!: b) = finite a && finite b
  bashZero (a :!: b) = bashZero a :!: bashZero b
  approxEq b (u :!: v) (x :!: y) = approxEq b u x && approxEq b v y

converge
  :: Approx t
  => Int -> Int -> [t] -> Maybe t
{-# SPECIALIZE converge :: Int -> Int -> [Complex Double] -> Maybe (Complex Double) #-}
{-# SPECIALIZE converge :: Int -> Int -> [Pair (Complex Double) (Complex Double)] -> Maybe (Pair (Complex Double) (Complex Double)) #-}
converge b n = lastMay . trim . filter finite . take n
  where
    trim [] = []
    trim [x] = [x]
    trim (x:ys@(y:_))
      | approxEq b x y = [x,y]
      | otherwise = trim ys

class Num t => Square t where
  sqr :: t -> t
  sqr x = x * x
  double :: t -> t
  double x = x + x

instance Square Double

instance (RealFloat t, Square t) => Square (Complex t) where
  {-# SPECIALIZE instance Square (Complex Double) #-}

magnitudeSquared :: Square r => Complex r -> r
{-# SPECIALIZE magnitudeSquared :: Complex Double -> Double #-}
magnitudeSquared (x :+ y) = sqr x + sqr y

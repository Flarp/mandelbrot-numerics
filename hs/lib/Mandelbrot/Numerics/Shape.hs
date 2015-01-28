{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Shape
  ( Shape(..)
  , shape
  ) where

import Mandelbrot.Numerics.Complex

data Shape = Cardioid | Circle
  deriving (Eq, Ord, Enum, Bounded, Read, Show)

shape
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Maybe Shape
shape p c
  | p < 1 || notFiniteC e = Nothing
  | magnitudeSquared e < magnitudeSquared (e - 1) = Just Cardioid
  | otherwise = Just Circle
  where
    e = go 1 c 1 1 0 0
    go i !z !dc !dz !dcdc !dcdz
      | i == p = - (dcdc / (double dc) + dcdz / dz) / (dc * dz)
      | otherwise = go
          (i + 1) (sqr z + c) (double (z * dc) + 1) (double (z * dz))
          (double (z * dcdc + sqr dc)) (double (z * dcdz + dc * dz))

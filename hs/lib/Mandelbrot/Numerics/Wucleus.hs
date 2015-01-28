{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Wucleus
  ( wucleus
  ) where

import Mandelbrot.Numerics.Complex

wucleus
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Complex r -> [Complex r]
wucleus p c0 z0
  | p < 1 || notFiniteC c0 || notFiniteC z0 = []
  | otherwise = go 0 z0 1
  where
    go i !z !dz
      | i == p = z' : wucleus p c0 z'
      | otherwise = go (i + 1) (sqr z + c0) (double (z * dz))
      where
        z' = z0 - (z - z0) / (dz - 1)

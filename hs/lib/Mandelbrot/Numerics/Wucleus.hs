{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Wucleus
  ( wucleus
  ) where

import Mandelbrot.Numerics.Complex

wucleus
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> Complex r -> [Complex r]
{-# SPECIALIZE wucleus :: Int -> Complex Double -> Complex Double -> [Complex Double] #-}
wucleus p c0 z0
  | p < 1 || notFinite c0 || notFinite z0 = []
  | otherwise = go 0 z0 1
  where
    go i !z !dz
      | i == p = z' : wucleus p c0 z'
      | otherwise = go (i + 1) (sqr z + c0) (double (z * dz))
      where
        z' = z0 - (z - z0) / (dz - 1)

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Nucleus
  ( nucleus
  ) where

import Mandelbrot.Numerics.Complex

nucleus
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> [Complex r]
{-# SPECIALIZE nucleus :: Int -> Complex Double -> [Complex Double] #-}
nucleus !p !c0@(cx :+ cy)
  | p > 0 && fin cx && fin cy = go 0 0 0
  | otherwise = []
  where
    fin x = not (isNaN x) && not (isInfinite x)
    go !i !z !dc
      | i == p = c' : nucleus p c'
      | otherwise = go (i + 1) (sqr z + c0) (double (z * dc) + 1)
      where
        c' = c0 - z / dc

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Misiurewicz
  ( misiurewicz
  , misiurewiczNaive
  ) where

import Data.List
  ( foldl'
  )

import Data.Strict.Tuple
  ( Pair((:!:))
  )

import Mandelbrot.Numerics.Complex

misiurewiczNaive
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Int -> Complex r -> [Complex r]
{-# SPECIALIZE misiurewiczNaive :: Int -> Int -> Complex Double -> [Complex Double] #-}
misiurewiczNaive pp p c0
  | pp > 0 && p > 0 && finite c0 = go 0 0 0 0 0
  | otherwise = []
  where
    go i !zp !dcp !z !dc
      | i == pp + p = c' : misiurewiczNaive pp p c'
      | i == pp = go (i + 1) z dc (sqr z + c0) (double (z * dc) + 1)
      | otherwise = go (i + 1) zp dcp (sqr z + c0) (double (z * dc) + 1)
      where
        c' = c0 - (z - zp) / (dc - dcp)

misiurewicz
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Int -> Complex r -> [Complex r]
{-# SPECIALIZE misiurewicz :: Int -> Int -> Complex Double -> [Complex Double] #-}
misiurewicz pp p c0
  | pp > 0 && p > 0 && finite c0 = go 0 [] 0 0
  | otherwise = []
  where
    go i !ps !z !dc
      | i == pp + p = c' : misiurewicz pp p c'
      | otherwise = go (i + 1) ((z:!:dc):ps) (sqr z + c0) (double (z * dc + 1))
      where
        c' = c0 - f / df
        f = g / h
        df = (dg * h - g * dh) / (h * h)
        g = z - zp
        dg = dc - dcp
        zp :!: dcp = ps !! (p - 1)
        h :!: dh0 = foldl' prodsum (1 :!: 0) $ zipWith sub ps (drop p ps)
        dh = h * dh0
        sub (xp :!: xs) (yp :!: ys) = dp :!: ((xs - ys) / dp)
          where dp = xp - yp
        prodsum (xp :!: xs) (yp :!: ys) = (xp * yp) :!: (xs + ys)

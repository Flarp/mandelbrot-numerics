{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.BoxPeriod
  ( boxPeriod
  ) where

import Mandelbrot.Numerics.Complex
import Mandelbrot.Numerics.Progress

cross
  :: Num r
  => Complex r -> Complex r -> r
cross (ax :+ ay) (bx :+ by) = ay * bx - ax * by

crossesPositiveRealAxis
  :: RealFloat r
  => Complex r -> Complex r -> Bool
crossesPositiveRealAxis a@(_:+ay) b@(_:+by)
  | signum ay /= signum by = s == t
  | otherwise = False
  where
    d@(_:+dy) = b - a
    s = signum dy
    t = signum (d `cross` a)

surroundsOrigin
  :: RealFloat r
  => Complex r -> Complex r -> Complex r -> Complex r -> Bool
surroundsOrigin a b c d
  = odd . length . filter id
  $ zipWith crossesPositiveRealAxis [a,b,c,d] [b,c,d,a]

boxPeriod
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Complex r -> r -> Progress Int
boxPeriod c r = go 1 c0 c1 c2 c3
  where
    r' = negate r
    c0 = c + (r' :+ r')
    c1 = c + (r  :+ r')
    c2 = c + (r  :+ r )
    c3 = c + (r' :+ r )
    go !p z0 z1 z2 z3
      | notFinite z0 = Failed p
      | notFinite z1 = Failed p
      | notFinite z2 = Failed p
      | notFinite z3 = Failed p
      | surroundsOrigin z0 z1 z2 z3 = Done p
      | otherwise = Continue p $
          go (p + 1) (sqr z0 + c0) (sqr z1 + c1) (sqr z2 + c2) (sqr z3 + c3)

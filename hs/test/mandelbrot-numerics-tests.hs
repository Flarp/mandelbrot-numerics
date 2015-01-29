{-# LANGUAGE TemplateHaskell #-}
module Main ( main ) where

import Data.Ratio
  ( (%)
  , denominator
  )
import Data.Word
  ( Word8
  , Word16
  )
import System.Exit
  ( exitFailure
  , exitSuccess
  )
import Test.QuickCheck
  ( Arbitrary(..)
  , arbitrarySizedBoundedIntegral
  , choose
  , vector
  , elements
  , Property
  , property
  , (==>)
  , quickCheckAll
  )

import Mandelbrot.Numerics

data Island = Island !Int !(Complex Double)
  deriving Show

islands :: [Island]
islands =
  [ Island 1 0
  , Island 3 (-1.754877666246693)
  , Island 4 ((-0.15652016683375508) :+ 1.0322471089228318)
  , Island 4 (-1.9407998065294847)
  , Island 4 ((-0.15652016683375508) :+ (-1.0322471089228318))
  ]

newtype IA8 = IA8{ unIA8 :: Rational } deriving Show
newtype IA16 = IA16 Rational deriving Show
newtype IA8s = IA8s [IA8] deriving Show

instance Arbitrary Island where
  arbitrary = elements islands

instance Arbitrary IA8 where
  arbitrary = do
    a <- arbitrarySizedBoundedIntegral
    b <- arbitrarySizedBoundedIntegral
    let n = min a b
        d = max a b
        fi :: Word8 -> Integer
        fi = fromIntegral
    return . IA8 $ (fi n + 1) % (fi d + 1)

instance Arbitrary IA16 where
  arbitrary = do
    a <- arbitrarySizedBoundedIntegral
    b <- arbitrarySizedBoundedIntegral
    let n = min a b
        d = max a b
        fi :: Word16 -> Integer
        fi = fromIntegral
    return . IA16 $ (fi n + 1) % (fi d + 1)

instance Arbitrary IA8s where
  arbitrary = do
    len <- choose (1, 10)
    IA8s `fmap` vector len

prop_parent_child :: Property
prop_parent_child = property $ \(IA16 s) (Island np n) ->
  0 < s && s < 1 ==>
  case child np n s of
    Just (Child p c) -> case parent p c of
      Just (Parent _ t _ _) -> s == t
      _ -> False
    _ -> False

prop_parents_children :: Property
prop_parents_children = property $ \(IA8s xs) (Island np n) ->
  let ss = map unIA8 xs
      cs = children np n ss
      Child p c = last cs
      ps = parents p c
      ts = reverse [ t | Parent _ t _ _ <- ps ]
  in  all (\s -> 0 < s && s < 1) ss &&
      product (map denominator ss) * fromIntegral np <= 2^(20 :: Int) ==>
      ss == ts

return []
main :: IO ()
main = do
  r <- $quickCheckAll
  if r then exitSuccess else exitFailure

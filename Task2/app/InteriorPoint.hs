{-# OPTIONS_GHC -Wno-orphans #-}
module InteriorPoint (interiorPointWithAccuracy, interiorPointIteration, dot) where

import Data.Matrix
import Data.Vector
import Data.Either
import Prelude hiding (sum, null, zipWith, length, map, minimum, foldl1, foldl, replicate)

inverse' :: (Fractional a, Eq a) => Matrix a -> Matrix a
inverse' = fromRight (error "matrix is not invertable") . inverse

instance Num a => Num (Vector a) where
  (*) :: Num a => Vector a -> Vector a -> Vector a
  (*) = undefined

  (+) :: Num a => Vector a -> Vector a -> Vector a
  (+) = zipWith (+)

  abs :: Num a => Vector a -> Vector a
  abs = map abs

  signum :: Num a => Vector a -> Vector a
  signum = map signum

  fromInteger :: Num a => Integer -> Vector a
  fromInteger = singleton . fromInteger

  negate :: Num a => Vector a -> Vector a
  negate = map negate

dot :: Num a => Vector a -> Vector a -> a
dot v = sum . zipWith (*) v

interiorPointIteration :: (Num a, Ord a, Fractional a, Eq a) =>
                           a -> Matrix a -> Vector a -> Vector a -> Vector a
interiorPointIteration alpha a c x =
  let d = diagonal 0 x
      a' = a * d
      c' = d * colVector c
      i = identity (ncols a)
      a't = transpose a'
      p = i - a't * inverse' (a' * a't) * a'
      c_p = getMatrixAsVector $ p * c'
      nu = minimum c_p
      x' = replicate (length c_p) 1 + map (* (alpha / abs nu)) c_p
      nextX = getMatrixAsVector $ d * colVector x'
  in nextX

interiorPointWithAccuracy :: (Num a, Ord a, Fractional a, Eq a, Floating a) =>
                              a -> a -> Matrix a -> Vector a -> Vector a -> Vector a
interiorPointWithAccuracy eps alpha a c x =
  let iterations = iterate (interiorPointIteration alpha a c) x
  in solutionWithAccuracy eps iterations

norm :: (Num a, Floating a, Fractional a) => Vector a -> a
norm v | null v = 0
norm v = sqrt $ foldl1 sumOfSquares v
  where sumOfSquares a b = a^^(2 :: Integer) + b^^(2 :: Integer)

solutionWithAccuracy :: (Ord a, Num a, Floating a, Fractional a) => a -> [Vector a] -> Vector a
solutionWithAccuracy eps (x1:x2:_) | norm (x1 - x2) <= eps = x1
solutionWithAccuracy eps (_:xs) = solutionWithAccuracy eps xs
solutionWithAccuracy _ [] = undefined


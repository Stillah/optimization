module Main where

import System.Process (createProcess, std_out, std_in, StdStream(CreatePipe), proc)
import Data.Matrix
import qualified Data.Matrix as Matrix
import qualified Data.Vector as Vector
import Data.Vector (Vector)
import GHC.IO.Handle (hClose, hGetContents)
import System.IO (hPrint)
import Data.Vector.Generic (foldl')
import Text.Printf (printf, PrintfArg)
import Data.List (dropWhileEnd)
import Text.Read (readMaybe)
import GHC.TopHandler (flushStdHandles)
import Control.Exception (handle)
import System.IO.Error (isDoesNotExistError)

import InteriorPoint

precisionOf :: (Show a, Real a) => a -> Int
precisionOf = length . dropWhile (/= '.') . show . abs

numberFormat :: (Show a, Real a, PrintfArg a) => a -> a -> String
numberFormat eps a =
  let precision = precisionOf eps
      formatted = printf ("%." <> show precision <> "f") a
  in if '.' `elem` formatted
    then dropWhileEnd (== '.') $ dropWhileEnd (== '0') formatted
    else formatted

outputFormat :: (Show a, PrintfArg a, Real a) => a -> Vector a -> a -> String
outputFormat eps solutionVec optimum =
  let prettySolutionVec = foldl' (\str a -> str<>" "<>numberFormat eps a) "" solutionVec
  in
  "Vector of decision variables x*:\n"<>
  prettySolutionVec <>
  "\nOptimum value of the objective function:\n "<>
  numberFormat eps optimum

readInt :: String -> IO Int
readInt prompt = do
  putStrLn prompt
  input <- getLine
  case readMaybe input of
    Just n -> return n
    _ -> putStrLn "Invalid number. Please enter an integer." >> readInt prompt

readVectorOfSize :: Read a => Int -> String -> IO (Vector a)
readVectorOfSize len prompt = do
  putStrLn prompt
  input <- getLine
  let parsed = mapM readMaybe (words input)
  case parsed of
    Just values | length values == len -> return (Vector.fromList values)
    _ -> do
      putStrLn "Invalid input. Please try again"
      readVectorOfSize len prompt

readVectorThat :: Read a => (Vector a -> Bool) -> String -> IO (Vector a)
readVectorThat p prompt = do
  putStrLn prompt
  input <- getLine
  let parsed = mapM readMaybe (words input)
  case parsed of
    Just values | p (Vector.fromList values) -> return (Vector.fromList values)
    _ -> do
      putStrLn "Invalid input. Please try again"
      readVectorThat p prompt

readMatrixOfSize :: Read a => Int -> Int -> String -> IO (Matrix a)
readMatrixOfSize m n prompt = do
  putStrLn prompt
  putStr $ show (m * n)
  putStr " values: "
  flushStdHandles
  input <- getLine
  case mapM readMaybe (words input) of
    Just vals | length vals == m * n -> return (Matrix.fromList m n vals)
    _ -> putStrLn "Invalid matrix input. Please try again" >> readMatrixOfSize m n prompt

readFractional :: (Read a, Fractional a) => String -> IO a
readFractional prompt = do
  putStrLn prompt
  input <- getLine
  case readMaybe input of
    Just eps -> return eps
    _ -> putStrLn "Invalid input. Please enter a valid number." >> readFractional prompt

readInputs :: (Read a, Fractional a, Eq a, Ord a) =>
               IO (Vector a, Matrix a, Vector a, Vector a, a)
readInputs = do
  m <- readInt "Enter number of constraint equations:"
  n <- readInt "Enter number of variables:"
  c <- readVectorOfSize n "Enter vector of coefficients for the objective function - c:"
  a <- readMatrixOfSize m n "Enter matrix of coefficients of constraint function - A"
  b <- readVectorOfSize m "Enter vector of right-hand side numbers - b:"
  x <- readVectorThat (\x -> length x == n
                          && all (> 0) x
                          && a * colVector x == colVector b)
    "Enter your initial starting point vector - x:"
  eps <- readFractional "Enter approximation accuracy epsilon:"
  return (c, a, x, b, eps)

printInteriorPointSolution :: (Ord a, Real a, Show a, PrintfArg a, Floating a) =>
                               a -> a -> Matrix a -> Vector a -> Vector a -> IO ()
printInteriorPointSolution eps alpha a c x = do
  let solutionVec = interiorPointWithAccuracy eps alpha a c x
  let optimum = c `dot` solutionVec

  putStrLn $ "\nAlpha = " <> show alpha <> ":"
  putStrLn $ outputFormat eps solutionVec optimum

printSimplexSolution :: (Read a, Ord a, Show a, Real a) =>
                         a -> Matrix a -> Vector a -> Vector a -> IO ()
printSimplexSolution eps a b c = handle handler $ do
  let n = ncols a
  let m = nrows a

  (Just stdin, Just stdout, _, _) <- createProcess (proc "./simplex" [])
    { std_in = CreatePipe, std_out = CreatePipe }

  let send :: Show a => a -> IO ()
      send = hPrint stdin
  send n
  send m
  mapM_ send c
  mapM_ send (a <|> colVector b)
  send eps

  hClose stdin

  [optimumStr, solutionVecStr] <- lines <$> hGetContents stdout
  let solutionVec = Vector.fromList $ map read $ words solutionVecStr
  let optimum = read optimumStr :: Double

  putStrLn "\nSimplex method:"
  putStrLn $ outputFormat (realToFrac eps) solutionVec optimum
  where
    handler :: IOError -> IO ()
    handler err
      | isDoesNotExistError err = putStrLn "Error: executable for simplex method does not exist"
      | otherwise = ioError err

main :: IO ()
main = do
  (c, a, x, b, eps) <- readInputs

  let alpha1 = 0.5 :: Double
  printInteriorPointSolution eps alpha1 a c x
  let alpha2 = 0.9 :: Double
  printInteriorPointSolution eps alpha2 a c x
  printSimplexSolution eps a b c


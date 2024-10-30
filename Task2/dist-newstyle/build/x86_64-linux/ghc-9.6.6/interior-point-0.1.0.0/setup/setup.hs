import Distribution.Simple
import Distribution.Simple.Setup
import Distribution.Types.HookedBuildInfo
import System.Process (callProcess)

main :: IO ()
main = defaultMainWithHooks simpleUserHooks {
  preBuild = \_ _ -> do
    callProcess "g++" ["-o", "simplex", "../Task1/Simplex.cpp", "-I../Task1"]
    return emptyHookedBuildInfo
}

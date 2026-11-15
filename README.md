# An improved lower bound for piercing the chessboard by lines

Ivan Matic and Rados Radoicic

This repository contains the verifier for the computer-assisted part of the paper above
(Theorem 3.2, statements (V1), (V2), (V3) of Section 4).

## Building

    make

On macOS with a Homebrew-installed GCC, use `make forMac` instead of `make`.

## Running

    ./piercing -i W14.txt

Invoked without arguments, the program prints its contract. The run takes about one second
and ends with the line `THE THEOREM WAS PROVED`.

## Third-party code

The directory `src/boost` contains an unmodified vendored subset of Boost 1.91.0
(Multiprecision and Config), used for fixed-width exact integer arithmetic.

## License

MIT for the authors' code; the vendored Boost subset in `src/boost` is under the Boost
Software License 1.0. See `LICENSE.txt`.

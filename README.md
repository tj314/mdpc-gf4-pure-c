# QC-MDPC McEliece cryptosystem over GF(4)

This project contains an implementation of [QC-MDPC McEliece cryptosystem over GF(4)](https://ieeexplore.ieee.org/document/8893339/). This is a research-only implementation and is not safe against side-channel attacks. DO NOT use this in real software.

This project was implement as a part of my master's thesis at FEI STU in Bratislava.

## How to run

You can either use an IDE such as Clion or you can run from terminal:

```bash
mkdir build
cd build
cmake ..
cmake --build .
./mdpc-gf4
```

Additionally, one may specify a build mode:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cmake --build .
./mdpc-gf4
```

The default mode is `Release`. There is also a `Testing` mode to run unit tests.

## What is implemented?

- finite field GF(4)

- polynomials and polynomial operations over GF(4)

- random generation of vectors

- key generating

- encoding and encryption

- decoding (with multiple decoders) and decryption

## Contributions

We will gladly accept your contributions! Feel free to create issues, forks, MRs... Try to use [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/) as much as possible.



























































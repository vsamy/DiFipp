DiFi++
======

[License BSD 2-Clause](https://tldrlegal.com/license/bsd-2-clause-license-(freebsd)#fulltext)

DiFi++ is a small c++ header-only library  for **DI**gital **FI**lters based on rational transfer functions such as the butterworth filter and the moving average filter. DiFi++ is using the Eigen library for math computations.

The implementation is based on well written article from Neil Robertson.
Please check out the followings

* [Butterworth filter](https://www.dsprelated.com/showarticle/1119.php)

* [Highpass filters](https://www.dsprelated.com/showarticle/1135.php)

* [Bandpass filters](https://www.dsprelated.com/showarticle/1128.php)

* [Band-reject filters](https://www.dsprelated.com/showarticle/1131.php)

and the differentiators from [Pavel Holoborodko](http://www.holoborodko.com/pavel/)

The library has been tested against Matlab results.

A doxygen documentation is generated when compiling.

Installing
-----

This is an header-only library so there is nothing to compile (but the documentation)

```bash
git clone --recursive https://github.com/vsamy/DiFi++
cd DiFi++
mkdir build
cd build
cmake ..
make install
```

Note
-----

The method used is close but somewhat different from Matlab methods and Butterworth band-reject has quite different results (precision of 1e-8).
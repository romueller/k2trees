# k2trees

Library of _k2-trees_ and variations of them.

_k2trees_ were proposed in [Compact representation of Web graphs with extended functionality](http://doi.org/10.1016/j.is.2013.08.003) (Brisaboa et al. 2014).

The implementation of the original _k2-tree_ (mostly) follows the description in above paper, but differs, e.g., in the incorporated optimisations.

The new variations generalise / transfer the general concept of _k2-trees_ to genuinely rectangular cases.

## Installation
Get the source code from GitHub and compile the library:

```sh
git clone https://github.com/romueller/k2trees.git
cd k2trees
make
```

To install the library into the `include` and `lib` subdirectories of `/usr/local`, run:

```sh
make install
```

A different prefix can be specified via `INSTALL_PREFIX=/path/to/directory`.


To remove the library again, run:

```sh
make uninstall
```


## Required software
 * C++ (GCC 4.9.2 or higher)
 * [SDSL](https://github.com/simongog/sdsl-lite)
 * make

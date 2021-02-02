# cohortid

Identifies cohorts based on the the presence or absence of some trait.

## Install

Automatic installation has only been tested on Ubuntu 18.04 and
macOS High Sierra, though it should work without issue on any Unix-alike.


Make sure you have [Julia](https://julialang.org/) in your `$PATH`.

```julia
pkg> add https://github.com/dorbarker/cohortid
```

This will install `cohortid` in `~/.local/bin`. If you have not done so already,
make sure this directory is in your `$PATH`.

```sh
echo "$HOME/.local/bin:$PATH" >> ~/.bashrc
source ~/.bashrc
```

## Usage

```sh
usage: cohortid --size-1 SIZE-1 --size-2 SIZE-2 --pos-1 POS-1
                --pos-2 POS-2 --neg-1 NEG-1 --neg-2 NEG-2
                --clusters CLUSTERS --metadata METADATA
                --variable VARIABLE --positive-value POSITIVE-VALUE
                -o OUTPUT [-h]

optional arguments:
  --size-1 SIZE-1       The lower bound size for consideration (type:
                        Int64)
  --size-2 SIZE-2       (type: Int64)
  --pos-1 POS-1         (type: Float64)
  --pos-2 POS-2         (type: Float64)
  --neg-1 NEG-1         (type: Float64)
  --neg-2 NEG-2         (type: Float64)
  --clusters CLUSTERS
  --metadata METADATA
  --variable VARIABLE
  --positive-value POSITIVE-VALUE

  -o, --output OUTPUT
  -h, --help            show this help message and exit
```

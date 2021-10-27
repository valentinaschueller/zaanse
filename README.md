# zaaNSE: A Held & Hou Atmosphere

zaaNSE is the attempt to simulate a zonally averaged atmosphere, using the Navier-Stokes equations -- it is also the failed attempt to come up with a good name for a toy project, which is why it might come up when someone is looking at [areas in the Netherlands](https://en.wikipedia.org/wiki/Zaanse_Schans).

This is an experiment to see if we can implement a stable 2D Finite Difference discretization of the zonally averaged 3D Navier-Stokes equations on the sphere.
Ideally, we get something which produces similar results as Held & Hou in their 1980 paper. \[1\]
And if not, we will for sure have learned something ;)

This is a continuation of the work we did at the [Ferienakademie 2021](http://www.ferienakademie.de/) and which I prepared theoretically in [this repository](https://github.com/valentinaschueller/ferienakademie-2021-presentation).

\[1\]: Held, Isaac M., and Arthur Y. Hou. “Nonlinear Axially Symmetric Circulations in a Nearly Inviscid Atmosphere.” Journal of the Atmospheric Sciences 37, no. 3 (March 1, 1980): 515–33. https://doi.org/10.1175/1520-0469(1980)037<0515:NASCIA>2.0.CO;2.

## Installing & Running zaaNSE

I recommend that you use a fresh Python environment, e.g. using Pipenv or Miniconda.
For ease of use, I also recommend that you install NumPy and Matplotlib before running `pip install`.

You can install zaaNSE best using 

```bash
$ git clone ...
$ cd zaanse
$ pip install -e .
```

You can run zaaNSE from the terminal using

```bash
$ zaanse -T 5 -dt 1 -v
$ # run zaanse -h for help and explanation of CLI parameters
```

To run the tests, you need to have [pytest](https://pytest.org) installed.
Run `pytest .` in the top level directory of the package.
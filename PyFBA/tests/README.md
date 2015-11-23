# Unit tests for PyFBA

The files in this directory are a collection of unit tests that test the accuracy of the code. If you have [installed
GLPK and PyGLPK (as we describe)](../INSTALLATION.md) then you should be able to test that installation with:

```
nosetest tests/test_lp.py
```

This will run a single test of the linear programming solver.

If you have installed PyFBA correctly, you should be able to run:

```
nosetest tests/
```

and test the entire code base.

Please let us know of any issues with the code.
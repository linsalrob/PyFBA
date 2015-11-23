# Linear Program Solvers

This package has wrappers for linear programming solvers. The essential methods that must be defined by the solver are:

* Load

```
    def load(matrix):
    or
    def load(matrix, rowheaders=None, colheaders=None):
    Load a matrix, optionally including the row and column headers. The matrix should be a two-dimensional matrix (a 
    list of lists)
```

* Row bounds:

```
    def row_bounds(bounds):
    Accept a list of tuples that define the row bounds
```

* Columns bounds:

```
    def col_bounds(bounds):
    Accept a list of tuples that define the column bounds
```

* Objective coefficient

```
    def objective_coefficients(coeff):
    A method that accepts a list that represents the objective coefficient of the problem. For FBA, this is usually the
    biomass equation.
```
    
* Solve

```
    def solve():
    Solve the problem and return the status of the solver and the value of the solution (i.e. the flux through the
    objective coefficient).
```

You can replace the solver used in PyFBA with a solver of your choice. We have used the GNU Linear Programming Toolkit
because it is freely available and compatible with all systems. However, it is not the fastest solver available, and so
you may prefer to replace it. We are working on wrappers for other solvers and will release them. In the meantime, if
you implement the above methods and change the first import line of [__init__.py](__init__.py) your wrapper should work.
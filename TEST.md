# Testing PyFBA

PyFBA comes with a complete test suite that tests different aspects of the code. We use [nosetests](http://nose.readthedocs.io/en/latest/) that you should be able to easily install on your machine. Once you have nosetests installed, you should be able to run:

```
cd PyFBA
nosetests tests
```

We have updated PyFBA to work with [Python 3](https://www.python.org/download/releases/3.0/), and you should also be able to test that using nosetests:

```
cd PyFBA
nosetests3 tests
```

If there are errors in the tests, please post an [issue on GitHub](https://github.com/linsalrob/PyFBA/issues) so we can address it.

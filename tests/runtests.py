#!/usr/bin/python

# Runner code for proteomeScoutAPI test suite
#
# To use, simply run
#
# python runtests.py from your terminal
# 

import unittest
import __init__ as test

# build and run the test suite
suite = unittest.TestLoader().loadTestsFromTestCase(test.unit_tests.TestFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)

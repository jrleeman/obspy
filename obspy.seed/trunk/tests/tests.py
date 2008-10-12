# -*- coding: utf-8 -*-

import doctest
import unittest

def suite():
    suite = unittest.TestSuite()
    suite.addTest(doctest.DocFileSuite('test_blockette010.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette011.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette030.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette033.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette034.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette050.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette052.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette053.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette054.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette057.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette058.txt'))
    suite.addTest(doctest.DocFileSuite('test_blockette061.txt'))
    return suite


if __name__ == '__main__':
    unittest.main(defaultTest='suite')
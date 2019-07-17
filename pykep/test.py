# -*- coding: utf-8 -*-

# Copyright 2017-2018 PaGMO development team
#
# This file is part of the PaGMO library.
#
# The PaGMO library is free software; you can redistribute it and/or modify
# it under the terms of either:
#
#   * the GNU Lesser General Public License as published by the Free
#     Software Foundation; either version 3 of the License, or (at your
#     option) any later version.
#
# or
#
#   * the GNU General Public License as published by the Free Software
#     Foundation; either version 3 of the License, or (at your option) any
#     later version.
#
# or both in parallel, as here.
#
# The PaGMO library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received copies of the GNU General Public License and the
# GNU Lesser General Public License along with the PaGMO library.  If not,
# see https://www.gnu.org/licenses/.

from __future__ import absolute_import as _ai

import unittest as _ut


class core_functions_test_case(_ut.TestCase):
    """Test case for the core functions

    """

    def runTest(self):
        self.run_epoch_test()

    def run_epoch_test(self):
        from .core import epoch
        epoch(julian_date=0, julian_date_type='mjd2000')


class lambert_test_case(_ut.TestCase):
    """Test case for the lambert problem class

    """

    def runTest(self):
        from .core import lambert_problem
        lambert_problem(r1=[1, 1, 0], r2=[0, 1, 0],
                        mu=1., cw=False, max_revs=0, tof=0.3)


class mga_1dsm_test_case(_ut.TestCase):
    """Test case for the mga1_dsm class

    """

    def runTest(self):
        self.run_construction_test()
        self.run_decode_times_and_vinf_test()

    def run_construction_test(self):
        from .trajopt import mga_1dsm
        # Correct use (nothrow)
        udp = mga_1dsm()
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='direct', tof=[[20, 400], [20, 400]])
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='eta', tof=500)
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='alpha', tof=[20, 500])
        self.assertEqual(udp.get_nobj(), 1)
        udp = mga_1dsm(tof_encoding='direct', tof=[
                       [20, 400], [20, 400]], multi_objective=True)
        self.assertEqual(udp.get_nobj(), 2)

        # Incorrect use (raise)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=34)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=[[400], [20, 400]])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='direct', tof=[20, 400])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='eta', tof=[20, 400])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='eta', tof=[[20, 400], [20, 400]])
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='alpha', tof=4)
        with self.assertRaises(TypeError):
            mga_1dsm(tof_encoding='alpha', tof=[[20, 400], [20, 400]])

    def run_decode_times_and_vinf_test(self):
        from .trajopt import mga_1dsm
        # Alpha
        udp = mga_1dsm(tof_encoding='alpha', tof=[20, 500])
        x = [10] + [0., 0., 1., 1., 0.5] + [0.4, 1.3, 0.5, 0.5] + [50]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 25.)
        self.assertAlmostEqual(retval[0][1], 25.)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)
        # Eta
        udp = mga_1dsm(tof_encoding='eta', tof=50)
        x = [10] + [0., 0., 1., 1., 0.5] + [0.4, 1.3, 0.5, 0.5]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 25.)
        self.assertAlmostEqual(retval[0][1], 12.5)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)
        # Direct
        udp = mga_1dsm(tof_encoding='direct', tof=[[10, 400], [10, 400]])
        x = [10] + [0., 0., 1., 1., 123] + [0.4, 1.3, 0.5, 321]
        retval = udp._decode_times_and_vinf(x)
        self.assertAlmostEqual(retval[0][0], 123)
        self.assertAlmostEqual(retval[0][1], 321)
        self.assertAlmostEqual(retval[1], 0.)
        self.assertAlmostEqual(retval[2], 0.)
        self.assertAlmostEqual(retval[3], 1.)


class gym_test_case(_ut.TestCase):
    """Test case for the gym

    """

    def runTest(self):
        self.run_rosetta_test()
        self.run_cassini2_test()


    def run_rosetta_test(self):
        from .trajopt import gym
        udp = gym.rosetta
        x = [1.53488329e+03, 4.56388378e-01, 9.51717655e-01, 4.18212047e+03,
             4.32159299e-01, 3.65256539e+02, 5.03363275e+00, 2.38949977e+00,
             4.55746823e-01, 7.09999954e+02, 1.79894273e+00, 1.05000003e+00,
             6.09083347e-01, 2.60816142e+02, 4.95158968e+00, 3.16049580e+00,
             6.89049263e-01, 7.29775762e+02, 4.30823655e+00, 1.10842692e+00,
             4.16075410e-01, 1.84999995e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], 1371.4992595125018)

    def run_cassini2_test(self):
        from .trajopt import gym
        udp = gym.cassini2
        x = [-7.75699976e+02,  9.15777367e-01,  4.06442043e-01,  3.21309562e+03,
        6.81118341e-01,  1.62660490e+02, -1.58051063e+00,  1.28479507e+00,
        4.72699902e-01,  4.24319550e+02,  4.30475919e+00,  1.15739933e+00,
        2.55718252e-01,  5.44489098e+01, -1.54332794e+00,  1.27160729e+00,
        9.00000000e-01,  5.88481599e+02,  4.76774269e+00,  7.00000000e+01,
        1.00000000e-02,  2.20000000e+03]
        self.assertAlmostEqual(udp.fitness(x)[0], 8511.423044703053)

def run_test_suite(level=0):
    """Run the full test suite.

    This function will raise an exception if at least one test fails.

    Args:
        level(``int``): the test level (higher values run longer tests)

    """

    retval = 0
    suite = _ut.TestLoader().loadTestsFromTestCase(core_functions_test_case)
    suite.addTest(lambert_test_case())
    suite.addTest(mga_1dsm_test_case())
    suite.addTest(gym_test_case())


    test_result = _ut.TextTestRunner(verbosity=2).run(suite)

    if len(test_result.failures) > 0 or len(test_result.errors) > 0:
        retval = 1
    if retval != 0:
        raise RuntimeError('One or more tests failed.')
#!/usr/bin/env python3

import unittest

from astro import *

class AstroTest(unittest.TestCase):

        def testSin(self):
                self.assertAlmostEqual(90, asin(sin(90)), "sin 90")
                self.assertAlmostEqual(0,  asin(sin(0)),  "sin 0")
                self.assertAlmostEqual(31, asin(sin(31)), "sin 31")

        def testCos(self):
                self.assertAlmostEqual(90, acos(cos(90)), "cos 90")
                self.assertAlmostEqual(0,  acos(cos(0)),  "cos 0")
                self.assertAlmostEqual(62, acos(cos(62)), "cos 62")

        def testTan(self):
                self.assertAlmostEqual(45, atan(tan(45)), "tan 45")
                self.assertAlmostEqual(0,  atan(tan(0)),  "tan 0")
                self.assertAlmostEqual(31, atan(tan(31)), "tan 31")

        def testJd(self):
                j = jd(1983, 1, 18, 7 + (12/60))
                self.assertAlmostEqual(2445352.8, j)

        def testSiderealTime(self):
                theta = sidereal_time(1982, 1, 1, 0)
                h, min, sec = deg2sex(theta)
                self.assertEqual(h, 6,      "mean sidereal time Greenwich h")
                self.assertEqual(min, 41,   "mean sidereal time Greenwich min")
                self.assertAlmostEqual(sec, 17.3, places=1,
                                       msg="mean sidereal time Greenwich sec")
                local = theta + 11.60833 / 15
                h, min, sec = deg2sex(local)
                self.assertEqual(h, 7,      "mean sidereal time Munich h")
                self.assertEqual(min, 27,   "mean sidereal time Munich min")
                self.assertAlmostEqual(sec, 43.3, places=1,
                                       msg="mean sidereal time Munich sec")

        def testGeoEqua2geoHori(self):
                t = (18 - 16) * 15.0
                ah = geoEqua2geoHori(t, 48, 20)
                self.assertAlmostEqual(ah[0], 51.3375, places=4, msg="azimuth")
                self.assertAlmostEqual(ah[1], 53.0068, places=4, msg="elevation")

        def testOrbit2helioEcl(self):
                lb = orbit2helioEcl(210, 30, 20)
                self.assertAlmostEqual(lb[0], 238.4812, places=4)
                self.assertAlmostEqual(lb[1], -9.8466, places=4)

        def testHelioEcl2geoEcl(self):
                v = helioEcl2geoEcl(150, 0, 1, 100, 1.3, 5)
                self.assertAlmostEqual(v[0], 1.4692, places=4, msg="beta")
                self.assertAlmostEqual(v[1], 90.0258, places=4, msg="lamda")
                self.assertAlmostEqual(v[2], 4.424226, places=6, msg="delta")

        def testSphe2cart(self):
                beta = 60; lam = 120; r = 5
                xyz = sphe2cart(beta, lam, r)
                self.assertAlmostEqual(xyz[0], -1.25, places=2, msg="x")
                self.assertAlmostEqual(xyz[1], 2.165064, places=6, msg="y")
                self.assertAlmostEqual(xyz[2], 4.330127, places=6, msg="z")
                s = cart2sphe(xyz[0], xyz[1], xyz[2])
                self.assertAlmostEqual(s[0], beta)
                self.assertAlmostEqual(s[1], lam)
                self.assertAlmostEqual(s[2], r)

        def testCart2sphe(self):
                x = -1.0; y = -2.5; z = 0.5
                s = cart2sphe(x, y, z)
                self.assertAlmostEqual(s[0], 10.5197, places=4, msg="beta")
                self.assertAlmostEqual(s[1], 248.1986, places=4, msg="lamda")
                self.assertAlmostEqual(s[2], 2.738613, places=6, msg="r")
                c = sphe2cart(s[0], s[1], s[2])
                self.assertAlmostEqual(c[0], x)
                self.assertAlmostEqual(c[1], y)
                self.assertAlmostEqual(c[2], z)

        def testGeoEcl2geoEqua(self):
               v = geoEcl2geoEqua(50, 290)
               self.assertAlmostEqual(v[0], 284.357, places=3, msg="alpha")
               self.assertAlmostEqual(v[1], 27.55, places=2, msg="delta")

        def testDeg2sec(self):
                d = 121.135
                dms = deg2sex(d)
                self.assertEqual(dms[0], 121)
                self.assertEqual(dms[1], 8)
                self.assertAlmostEqual(dms[2], 6.0)
                self.assertAlmostEqual(sex2deg(dms[0], dms[1], dms[2]), d)

        def testPositionSun(self):
                ra, dec = calcPositionSun(jd(2020, 1, 1, 12))
                self.assertAlmostEqual(ra, 281.44595, places=4)
                self.assertAlmostEqual(dec, -23.0221, places=4)

        def testPositionMoon(self):
                ra, dec = calcPositionMoon(jd(2020, 1, 1, 18))
                self.assertAlmostEqual(ra, 357.511, places=1)
                self.assertAlmostEqual(dec, -6.686, places=1)


if __name__ == '__main__':
        unittest.main()

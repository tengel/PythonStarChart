#!/usr/bin/env python3

import unittest
import math
import kepler
import ephem
import datetime

class TestPlanets(unittest.TestCase):

    def compareDegrees(self, a, b, delta, msg):
        if abs(a - b) > 300:
            if a > 300:
                a -= 360
            if b > 300:
                b -=  360
        self.assertAlmostEqual(a, b, delta=delta, msg=msg)

    def doChecks(self, ref, dut, dateStr):
        if dut.name == "Earth":
            sunDistanceRef = ref.earth_distance
        else:
            sunDistanceRef = ref.sun_distance

        self.assertAlmostEqual(sunDistanceRef, dut.distance_sun,
                               delta=0.001,
                               msg="Distance Sun, %s, %s" % (dut.name, dateStr))
        self.assertAlmostEqual(math.degrees(ref.hlat), dut.helio_lat,
                               delta=0.1,
                               msg="Heliocentric latitude, %s, %s" % (dut.name, dateStr))
        self.compareDegrees(math.degrees(ref.hlon), dut.helio_lon, 0.8,
                            "Heliocentric longutide, %s, %s" % (dut.name, dateStr))
        if dut.name == "Earth":
            return

        self.assertAlmostEqual(ref.earth_distance, dut.distance_earth,
                               delta=0.1,
                               msg="Distance earth, %s, %s" % (dut.name, dateStr))
        self.compareDegrees(math.degrees(ref.a_ra), dut.ra, 1.8,
                            "Geocentric, Equatorial, RA, %s, %s" % (dut.name, dateStr))
        self.compareDegrees(math.degrees(ref.a_dec), dut.declination, 0.7,
                            "Geocentric, Equatorial, Dec, %s, %s" % (dut.name, dateStr))

    def runTests(self, ref, dut):
        for y in range(1980, 2051):
            for m in range(1, 13):
                date = datetime.datetime(y, m, 1)
                myEarth = kepler.Earth()
                myEarth.calcHeliocentric(date)

                ref.compute("%s/%s/01" % (y, m))
                dut.calcHeliocentric(date)
                if dut.name != "Earth":
                    dut.calcGeocentric(myEarth)
                dateStr = "%s/%s" % (y, m)
                self.doChecks(ref, dut, dateStr)



    def test_mercury(self):
        self.runTests(ephem.Mercury(), kepler.Mercury())

    def test_venus(self):
        self.runTests(ephem.Venus(), kepler.Venus())

    def test_earth(self):
        self.runTests(ephem.Sun(), kepler.Earth())

    def test_mars(self):
        self.runTests(ephem.Mars(), kepler.Mars())

    def test_jupiter(self):
        self.runTests(ephem.Jupiter(), kepler.Jupiter())

    def test_saturn(self):
        self.runTests(ephem.Saturn(), kepler.Saturn())

    def test_uranus(self):
        self.runTests(ephem.Uranus(), kepler.Uranus())

    def test_neptune(self):
        self.runTests(ephem.Neptune(), kepler.Neptune())

if __name__ == '__main__':
    unittest.main()

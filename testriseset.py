#!/usr/bin/env python3

import unittest
import astro
import kepler
import ephem
import datetime
import math

def timeDiff(td1, td2):
    diff =  min(abs((td1 - td2).total_seconds()),
                abs((td2 - td1).total_seconds()))
    if diff > datetime.timedelta(hours=23).total_seconds():
        diff -= datetime.timedelta(hours=24).total_seconds();
    return round(diff)

class TestRiseSet(unittest.TestCase):

    def do_tests(self, elevation, ephemObj, testFunction, deltaMinutes):
        date    = datetime.date(2010, 1, 1)
        endDate = datetime.date(2040, 12, 31)
        longitude = 9.49
        latitude = 51.31
        observer = ephem.Observer()
        observer.pressure = 0
        observer.horizon = math.radians(elevation)
        observer.lat = math.radians(latitude)
        observer.lon = math.radians(longitude)
        dtFmt = "%Y-%m-%d %H:%M:%S"
        while date < endDate:
            observer.date = date
            refRise = observer.next_rising(ephemObj).datetime()
            refSet = observer.next_setting(ephemObj).datetime()
            myRise, mySet = testFunction(longitude, latitude, date)
            self.assertLess(
                timeDiff(refRise, myRise), (deltaMinutes*60),
                msg="rise time '%s' at '%s': expected: '%s' is: '%s'" %
                (ephemObj.name, date, refRise.strftime(dtFmt),
                 myRise.strftime(dtFmt)))
            self.assertLess(
                timeDiff(refSet, mySet), (deltaMinutes*60),
                msg="set time '%s' at '%s': expected: '%s' is: '%s'" %
                (ephemObj.name, date, refSet.strftime(dtFmt),
                 mySet.strftime(dtFmt)))
            date += datetime.timedelta(days=1)


    def test_stars(self):
        elevation = -0.566667
        ephemObject = ephem.star("Rigel")
        testFunction = astro.calcRiseSet_star
        objRa = math.degrees(ephemObject._ra) / 15
        objDec = math.degrees(ephemObject._dec)
        def doCalc(longitude, latitude, date):
            return astro.calcRiseSet_star(longitude, latitude, date, objRa,
                                          objDec)
        self.do_tests(elevation, ephemObject, doCalc, 3)


    def test_sun(self):
        elevation = -0.83333
        ephemObject = ephem.Sun()
        testFunction = astro.calcRiseSet_sun
        self.do_tests(elevation, ephemObject, testFunction, 3)


    @unittest.skip
    def test_moon(self):
        elevation = 0.133333
        ephemObject = ephem.Moon()
        testFunction = astro.calcRiseSet_moon
        self.do_tests(elevation, ephemObject, testFunction, 10)


    def test_moon2(self):
        date = datetime.date(2010, 1, 2)
        longitude = 15.0
        latitude = 50.0
        elevation = 0.133333
        observer = ephem.Observer()
        observer.pressure = 0
        observer.horizon = math.radians(elevation)
        observer.lat = math.radians(latitude)
        observer.lon = math.radians(longitude)
        observer.date = date
        refRise = observer.previous_rising(ephem.Moon()).datetime()
        refSet = observer.next_setting(ephem.Moon()).datetime()
        def objPosFunc(t):
            ra, dec = astro.calcPositionMoon(astro.julian_date(t))
            return ra / 15, dec
        myRise, mySet = astro.calcRiseSet_moon(longitude, latitude, date)
        dtFmt = "%Y-%m-%d %H:%M:%S"
        #print("rise ref:", refRise.strftime(dtFmt))
        #print("     my: ", myRise.strftime(dtFmt))
        #print("set  ref:", refSet.strftime(dtFmt))
        #print("     my: ", mySet.strftime(dtFmt))
        self.assertLess(
            timeDiff(refRise, myRise), (7*60),
            msg="rise time 'Moon' at '%s': expected: '%s' is: '%s'" %
            (date, refRise.strftime(dtFmt), myRise.strftime(dtFmt)))
        self.assertLess(
            timeDiff(refSet, mySet), (7*60),
            msg="set time 'Moon' at '%s': expected: '%s' is: '%s'" %
            (date, refSet.strftime(dtFmt), mySet.strftime(dtFmt)))


    def test_planets(self):
        elevation = -0.566667
        def doCalc(longitude, latitude, date):
            return kepler.calcRiseSet_planet(longitude, latitude, date, planet)
        ephemObject = ephem.Mercury()
        planet = kepler.Mercury()
        self.do_tests(elevation, ephemObject, doCalc, 1)
        ephemObject = ephem.Venus()
        planet = kepler.Venus()
        self.do_tests(elevation, ephemObject, doCalc, 1)
        ephemObject = ephem.Mars()
        planet = kepler.Mars()
        self.do_tests(elevation, ephemObject, doCalc, 11)
        ephemObject = ephem.Jupiter()
        planet = kepler.Jupiter()
        self.do_tests(elevation, ephemObject, doCalc, 4)
        ephemObject = ephem.Saturn()
        planet = kepler.Saturn()
        self.do_tests(elevation, ephemObject, doCalc, 4)
        ephemObject = ephem.Uranus()
        planet = kepler.Uranus()
        self.do_tests(elevation, ephemObject, doCalc, 4)
        ephemObject = ephem.Neptune()
        planet = kepler.Neptune()
        self.do_tests(elevation, ephemObject, doCalc, 4)


if __name__ == '__main__':
    unittest.main()

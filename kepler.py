import math
import sys
import copy
import pathlib
from astro import *

class OrbitalElements:
    def __init__(self):
        self.MA = 0.0 # Mean anomaly, M (degrees)
        self.EC = 0.0 # Eccentricity, e
        self.IN = 0.0 # Inclination w.r.t XY-plane, i (degrees)
        self.OM = 0.0 # Longitude of Ascending Node, Omega (degrees)
        self.W  = 0.0 # Argument of Perifocus, w (degrees)
        self.om = 0.0 # longitude of perihelion, omega_bar = (OM + W) % 360
        self.A  = 0.0 # Semi-major axis, a (au)
        self.N  = 0.0 # Mean motion, n (degrees/day)


class Planet:
    def __init__(self, name):
        self.name = name
        self.oe = OrbitalElements()
        # heliocentric, ecliptical (date/time required)
        self.helio_lon = 0
        self.helio_lat = 0
        self.distance_sun = 0
        # geocentric, ecliptical (position of earth required)
        self.ecliptic_lat = 0
        self.ecliptic_lon = 0
        self.distance_earth = 0
        # geocentric, equatorial
        self.declination = 0
        self.ra = 0
        # horizontal (position of observer required)
        self.altitude = 0
        self.azimuth = 0

    def _calcHeliocentric(self):
        """
        Calculate the heliocentric, ecliptical coordinates.
        Sets helio_lon, helio_lat and distance_sun of the object.
        """
        E = calcEccentricAnomaly(self.oe.MA, self.oe.EC)
        true_anomaly  = calcTrueAnomaly(self.oe.EC, E)
        self.distance_sun = calcDistanceSun(self.oe.A, self.oe.EC, E)
        omega = self.oe.om - self.oe.OM
        u = omega + true_anomaly
        self.helio_lon, self.helio_lat = orbit2helioEcl(u, self.oe.OM, self.oe.IN)

    def calcGeocentric(self, earth):
        """
        Calculate the geocentric, ecliptical coordinates ecliptic_lat,
        ecliptic_lon and the geocentric, equatorial coordinates
        declination and ra.
        """
        beta, lam, delta = helioEcl2geoEcl(earth.helio_lon, earth.helio_lat,
                                           earth.distance_sun,
                                           self.helio_lon, self.helio_lat,
                                           self.distance_sun)
        self.ecliptic_lat = beta
        self.ecliptic_lon = lam
        self.distance_earth = delta
        alpha, delta = geoEcl2geoEqua(beta, lam)
        self.declination = delta
        self.ra = alpha


class PlanetCsv(Planet):
    def __init__(self, name, csvFile):
        super().__init__(name)
        self.elements = self._readElements(csvFile)

    def _readElements(self, oeFile):
        """
        Read orbital elements from CSV file oeFile and return a list of OrbitalElement objects.
        """
        oeList = []
        csv = open(pathlib.Path(__file__).parent / oeFile)
        for line in csv:
            if line.startswith("$$SOE"):
                break
        for line in csv:
            if line.startswith("$$EOE"):
                break
            parts = line.split(",")
            if len(parts) != 15:
                raise Exception("invalid length: %s" % len(parts))
            oe = OrbitalElements()
            oe.MA = float(parts[9])
            oe.EC = float(parts[2])
            oe.IN = float(parts[4])
            oe.OM = float(parts[5])
            oe.W  = float(parts[6])
            oe.om = (oe.OM + oe.W) % 360
            oe.A  = float(parts[11])
            oe.N  = float(parts[8])
            oe.JD = float(parts[0])
            oeList.append(oe)
        csv.close()
        return oeList

    def _findClosestOe(self, date):
        """
        Find closest matching orbital elements for date.
        Returns OrbitalElement object, and time difference in days.
        """
        last_diff = sys.float_info.max
        e = 0
        for e in self.elements:
            diff =  julian_date(date) - e.JD
            if abs(last_diff) < abs(diff):
                e = last_e
                diff = last_diff
                break
            last_e = e
            last_diff = diff
        return copy.copy(e), diff

    def calcHeliocentric(self, date):
        self.oe, dayDiff = self._findClosestOe(date)
        self.oe.MA += (self.oe.N * dayDiff)
        self._calcHeliocentric()


class Mercury(Planet):
    def __init__(self):
        super().__init__("Mercury")

    def calcHeliocentric(self, date):
        self.oe.MA = 174.7947 + 149472.5153  * julian_century(date)
        self.oe.EC = 0.205634 + 0.000020 * julian_century(date)
        self.oe.IN = 7.0048 + 0.0019 * julian_century(date)
        self.oe.OM = 48.331 + 1.185 * julian_century(date)
        self.oe.om = 77.4552 + 1.5555 * julian_century(date)
        self.oe.A  = 0.387099
        self._calcHeliocentric()

class Venus(Planet):
    def __init__(self):
        super().__init__("Venus")

    def calcHeliocentric(self, date):
        self.oe.MA = 50.4071 + 58517.8039 * julian_century(date)
        self.oe.EC = 0.006773 - 0.000048 * julian_century(date)
        self.oe.IN = 3.3946 + 0.0010 * julian_century(date)
        self.oe.OM = 76.680 + 0.900 * julian_century(date)
        self.oe.om = 131.5718 + 1.4080 * julian_century(date)
        self.oe.A  = 0.723332
        self._calcHeliocentric()

class Earth(Planet):
    def __init__(self):
        super().__init__("Earth")

    def calcHeliocentric(self, date):
        self.oe.MA = 357.5256 + 35999.0498 * julian_century(date)
        self.oe.EC = 0.016709 - 0.000042 * julian_century(date)
        self.oe.IN = 0.0
        self.oe.OM = 0.0
        self.oe.om = 102.9400 + 1.7192 * julian_century(date)
        self.oe.A  = 1.0
        self._calcHeliocentric()

class Mars(Planet):
    def __init__(self):
        super().__init__("Mars")

    def calcHeliocentric(self, date):
        self.oe.MA = 19.3879 + 19139.8585 * julian_century(date)
        self.oe.EC = 0.093405 - 0.000092 * julian_century(date)
        self.oe.IN = 1.8496 - 0.0007 * julian_century(date)
        self.oe.OM = 49.557 + 0.771 * julian_century(date)
        self.oe.om = 336.0590 + 0.4438 * julian_century(date)
        self.oe.A  = 1.523692
        self._calcHeliocentric()

class Jupiter(PlanetCsv):
    def __init__(self):
        super().__init__("Jupiter", "data/horizons_jupiter.csv")

class Saturn(PlanetCsv):
    def __init__(self):
        super().__init__("Saturn", "data/horizons_saturn.csv")

class Uranus(PlanetCsv):
    def __init__(self):
        super().__init__("Uranus", "data/horizons_uranus.csv")

class Neptune(PlanetCsv):
    def __init__(self):
        super().__init__("Neptune", "data/horizons_neptune.csv")


def calcEccentricAnomaly(meanAnomaly, eccentricity):
    """
    Calculates the eccentric anomaly in degrees from mean anomaly and
    eccentricity. Calculation is done iteratively until the precision of
    epsilon is reached.
    """
    epsilon = 0.00001
    E = meanAnomaly
    while True:
        Enext = (E -
                 ((meanAnomaly - E + (180/math.pi) * eccentricity * sin(E)) /
                  (eccentricity * cos(E) - 1)))
        if abs(E - Enext) < epsilon:
            return Enext
        E = Enext

def calcTrueAnomaly(eccentricity, eccentricAnomaly):
    """
    Calculate the true anomaly.
    """
    t = (math.sqrt((1 + eccentricity) / (1 - eccentricity)) *
         tan(eccentricAnomaly / 2))
    v = atan(t) * 2
    if v < 0:
        return 360 + v
    else:
        return v

def calcDistanceSun(semiMajorAxis, eccentricity, eccentricAnomaly):
    """
    Calculate the distance from the planet to the sun in AU.
    """
    return semiMajorAxis * (1 - eccentricity * cos(eccentricAnomaly))


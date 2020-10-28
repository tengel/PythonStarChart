"""
The astro module contains basic astronomical calculations of angles, time and
coordinate systems.

If not otherwise stated, the formulas are based on:
Montenbruck O.; Grundlagen der Ephemeridenrechnung; Spektrum Akademischer
Verlag, München, 7. Auflage (2005)


----
"""

import math
import datetime

def sin(x):
    """
    Return the sine of x (measured in degrees).
    """
    return math.sin(math.radians(x))

def cos(x):
    """
    Return the cosine of x (measured in degrees).
    """
    return math.cos(math.radians(x))

def tan(x):
    """
    Return the tangent of x (measured in degrees).
    """
    return math.tan(math.radians(x))

def asin(x):
    """
    Return the arc sine (measured in degrees) of x.
    """
    return math.degrees(math.asin(x))

def acos(x):
    """
    Return the arc cosine (measured in degrees) of x.
    """
    return math.degrees(math.acos(x))

def atan(x):
    """
    Return the arc tangent (measured in degrees) of x.
    """
    return math.degrees(math.atan(x))

def jd(sY, sM, sD, sUT=0):
    """
    Return the julian date. (valid after 15 Oct 1582)

    :param int sY: Year
    :param int sM: Month
    :param int sD: Day
    :param float sUT: Hours in Universal Time
    :rtype: float
    """
    if sM <= 2:
        y = sY -1
        m = sM + 12
    else:
        y = sY
        m = sM
    b = math.floor(y/400.0) - math.floor(y / 100.0)
    jd = math.floor(365.25 * y) + math.floor(30.6001 * (m+1)) + b + 1720996.5 + sD + sUT/24.0
    return jd

def julian_date(dt):
    """
    Return the julian date (as float) from a python datetime object.
    """
    return jd(dt.year, dt.month, dt.day,
              (dt.hour + dt.minute / 60 + dt.second / 3600))

def julian_century(dt):
    """
    Return centuries since J2000.0 (as float) from a python datetime object.
    """
    return (julian_date(dt) - 2451545.0) / 36525

def sidereal_time(sY, sM, sD, sUT):
    """
    Return the mean sidereal time of Greenwich in hours.

    :param int sY: Year
    :param int sM: Month
    :param int sD: Day
    :param float sUT: Hours in Universal Time
    :rtype: float
    """
    jd0 = jd(sY, sM, sD)
    theta = 6.664520 + 0.0657098244 * (jd0 - 2451544.5) + 1.0027379093 * sUT
    return theta % 24

def geoEqua2geoHori(t, phi, delta):
    """
    Convert geocentric, equatorial coordinates (right ascension, declination) to
    horizontal coordinates (azimuth, elevation).

    :param float t: Hour angle in degree.
                    (t=local_sidereal_time - right_ascension)
    :param float phi: Geographical latitude of observer. (degree)
    :param float delta: Declination (degree)
    :return: Horizontal coordinates (azimuth, elevation) in degree
    :rtype: list(float)
    """
    x = sin(phi) * cos(delta) * cos(t) - cos(phi) * sin(delta)
    y = cos(delta) * sin(t)
    z = sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(t)
    r = math.sqrt(x**2 + y**2 + z**2)
    p = math.sqrt(x**2 + y**2)
    beta = atan(z/p)
    lamda_helper = 2 * atan(y / ( abs(x) + math.sqrt(x**2 + y**2)))
    if x == 0 and y == 0:
        lamda = 0
    elif x >= 0 and y >= 0:
        lamda = lamda_helper
    elif x >= 0 and y < 0:
        lamda = 360 + lamda_helper
    elif x < 0:
        lamda = 180 - lamda_helper
    return lamda, beta

def orbit2helioEcl(u, Omega, i):
    """
    Convert coordinates of the orbital plane to heliocentric, ecliptic
    coordinates (heliocentric ecliptic latitude, heliocentric ecliptic
    longitude).

    :param float u: Angular position of the planet. Angle between ascending node
                    and planet. (u = argument_of_periapsis + true_anomaly)
    :param float Omega: Longitude of the ascending node (degrees).
    :param float i: Inclination. Angle between the orbital plane and the ecliptic.
    :return: Heliocentric ecliptic longitude (0..360),
             heliocentric ecliptic latitude (-90..+90)
    """
    b = asin(sin(u) * sin(i))
    l = acos(cos(u) / cos(b)) + Omega
    lo1 = acos(cos(u) / cos(b))
    lo2 = asin((sin(u)*cos(i))/cos(b))
    if lo2 < 0:
        l = 360.0 - (lo1 - Omega)
    if l > 360:
        l = l - 360
    return l, b

def helioEcl2geoEcl(eLon, eLat, eDist, pLon, pLat, pDist):
    """
    Convert heliocentric, ecliptic coordinates to geocentric, ecliptic coordinates.

    :param float eLon: Heliocentric, ecliptic longitude of earth (degree)
    :param float eLat: Heliocentric, ecliptic latitude of earth (degree)
    :param float eDist: Earth distance from sun (AU)
    :param float pLon: Heliocentric, ecliptic longitude of planet (degree)
    :param float pLat: Heliocentric, ecliptic latitude of planet (degree)
    :param float pDist: Distance from sun (AU)
    :return: Geocentric ecliptic coordinates
             lamda (ecliptic lon), beta (ecliptic lat), Delta (earth distance)
    :rtype: list(float)
    """
    x = pDist * cos(pLat) * cos(pLon) - eDist * cos(eLat) * cos(eLon)
    y = pDist * cos(pLat) * sin(pLon)   - eDist * cos(eLat) * sin(eLon)
    z = pDist * sin(pLat)                   - eDist * sin(eLat)
    return cart2sphe(x, y, z)

def sphe2cart(beta, lamda, r):
    """
    Convert spherical coordinates (beta, lambda, r) to cartesian coordinates
    (x, y, z).

    :return: Cartesian coordinates x, y, z.
    """
    x = r * cos(beta) * cos(lamda)
    y = r * cos(beta) * sin(lamda)
    z = r * sin(beta)
    return x, y, z

def cart2sphe(x, y, z):
    """
    Convert cartesian coordinates (x, y, z) to  spherical coordinate
    (beta, lambda, r).

    :return: Spherical coordinates beta, lambda, r
    """
    r = math.sqrt(x**2 + y**2 + z**2)
    p = math.sqrt(x**2 + y**2)
    if p == 0:
        if   z > 0:
            beta = 90
        elif z == 0:
            beta = 0
        elif z < 0:
            beta = -90
    else:
        beta = atan(z / p)
    phi = 2 * atan(y / (abs(x) + math.sqrt(x**2 + y**2)))
    if   x == 0 and y == 0:
        lamda = 0
    elif x >= 0 and y >= 0:
        lamda = phi
    elif x >= 0 and y < 0:
        lamda = 360 + phi
    elif x < 0:
        lamda = 180 - phi
    return beta, lamda, r

def geoEcl2geoEqua(beta, lamb):
    """
    Convert geocentric, ecliptic coordinates (ecl. lat, ecl. long) to geocentric
    equatorial coordinates (right ascension, declination).

    The inclination of the ecliptic is fixed for J2000.

    :param float beta: Geocentric ecliptic latitude beta (degree)
    :param float lamb: Geocentric ecliptic longitude lambda (degree)
    :return: Equatorial coordinates
             (right ascension alpha, declination delta) in degree
    :rtype: list(float)
    """
    epsilon = 23.4392916667
    x = cos(beta) * cos(lamb)
    y = cos(epsilon) * cos(beta) * sin(lamb) - sin(epsilon) * sin(beta)
    z = sin(epsilon) * cos(beta) * sin(lamb) + cos(epsilon) * sin(beta)
    polar = cart2sphe(x, y, z)
    return polar[1], polar[0]

def calcPositionSun(jd):
    """
    Calculate the geocentric, equatorial position (ra, dec) of the sun for a
    julian date.
    Source: https://en.wikipedia.org/wiki/Position_of_the_Sun

    :param float jd: Julian date
    :return: Geocentric, equatorial coordinates
             (right ascension alpha, declination delta) in degree
    """
    n = jd - 2451545.0
    L = (280.460 + 0.9856474 * n) % 360
    g = (357.528 + 0.9856003 * n) % 360
    lamb = L + 1.915 * sin(g) + 0.020 * sin(2 * g) # geoc., ecliptic lon
    return geoEcl2geoEqua(0, lamb)   # ra, dec

def calcPositionMoon(jd):
    """
    Calculate the geocentric, equatorial position (ra, dec) of the moon for a
    julian date.

    :param float jd: Julian date
    :return: Geocentric, equatorial coordinates
             (right ascension alpha, declination delta) in degree
    """
    T  = (jd - 2451545.0) / 36525
    L0 = (218.31665 + 481267.88134 * T - 0.001327 * T**2) % 360.0
    l  = (134.96341 + 477198.86763 * T + 0.008997 * T**2) % 360.0
    l_ = (357.52911 + 35999.05029 * T + 0.000154 * T**2) % 360.0
    F  = (93.27210 + 483202.01753 * T - 0.003403 * T**2) % 360.0
    D  = (297.85020 + 445267.11152 * T - 0.001630 * T**2) % 360.0
    L1 = (22640 * sin(l) + 769 * sin(2 * l) + 36 * sin(3 * l)
          -4586 * sin(l - 2 * D)
          +2370 * sin(2 * D)
          -668  * sin(l_)
          -412  * sin(2 * F)
          -212  * sin(2 * l - 2 * D)
          -206  * sin(l + l_ - 2 * D)
          +192  * sin(l + 2 * D)
          -165  * sin(l_ - 2 * D)
          +148  * sin(l - l_)
          -125  * sin(D)
          -110  * sin(l + l_)
          -55   * sin(2 * F - 2 * D))
    lamb = L0 + (L1 / 60 / 60) # geocentric, ecliptic longitude
    B = (18520 * sin(F + lamb - L0 + 0.114 * sin(2 * F) + 0.150 * sin(l_))
         -526  * sin(+F - 2 * D)
         +44   * sin(+l + F - 2 * D)
         -31   * sin(-l + F - 2 * D)
         -25   * sin(-2 * l + F)
         -23   * sin(+l_ + F - 2 * D)
         +21   * sin(-l + F)
         +11   * sin(-l_ + F - 2 * D))
    beta = B / 60 / 60 # geocentric, ecliptic latitude
    return geoEcl2geoEqua(beta, lamb)   # ra, dec

def deg2sex(degree):
    """
    Convert degree to sexagesimal deg, min, sec.
    """
    deg = int(degree)
    min = (degree - deg) * 60
    sec = (min - int(min)) * 60
    return (deg, int(min), sec)

def sex2deg(d, m, s):
    """
    Convert sexagesimal deg, min, sec to degree float.
    """
    return d + (m / 60) + (s / 60 / 60)

def calcRiseSet_star(longitude, latitude, date, objRa, objDec):
    """
    Calculate the rise time and set time of a star.

    :param float longitude: Geographical longitude of observer (degree).
    :param float latitude: Geographical latitude of observer (degree).
    :param datetime date: Date (python date object).
    :param float objRa: Right ascension of star (h).
    :param float objDec: Declination of star (degree).
    :return: Returns (rise time, set time) as python datetime objects.
    """
    objType = "star"
    elevation = -0.566667
    r = calcRiseSet(objType, elevation, longitude, latitude,
                    lambda t: (objRa, objDec),
                    date, 12, True, 2)
    s = calcRiseSet(objType, elevation, longitude, latitude,
                    lambda t: (objRa, objDec),
                    date, 12, False, 2)
    return r, s

def calcRiseSet_sun(longitude, latitude, date):
    """
    Calculate the rise time and set time of the Sun.

    :param float longitude: Geographical longitude of observer (degree).
    :param float latitude: Geographical latitude of observer (degree).
    :param datetime date: Date (python date object).
    :return: Returns (rise time, set time) as python datetime objects.
    """
    objType = "sun"
    elevation = -0.83333
    def objPosFunc(t):
        ra, dec = calcPositionSun(julian_date(t))
        return ra / 15, dec
    r = calcRiseSet(objType, elevation, longitude, latitude, objPosFunc,
                    date, 6, True, 2)
    s = calcRiseSet(objType, elevation, longitude, latitude, objPosFunc,
                    date, 18, False, 2)
    return r, s

def calcRiseSet_moon(longitude, latitude, date):
    """
    Calculate the rise time and set time of the Moon.

    :param float longitude: Geographical longitude of observer (degree).
    :param float latitude: Geographical latitude of observer (degree).
    :param datetime date: Date (python date object).
    :return: Returns (rise time, set time) as python datetime objects.
    """
    objType = "moon"
    elevation = 0.133333
    def objPosFunc(t):
        ra, dec = calcPositionMoon(julian_date(t))
        return ra / 15, dec
    r = calcRiseSet(objType, elevation, longitude, latitude, objPosFunc,
                    date, 12, True, 5)
    s = calcRiseSet(objType, elevation, longitude, latitude, objPosFunc,
                    date, 12, False, 5)
    return r, s

def calcRiseSet(objType, elevation, longitude, latitude, objPosFunc, date, T0,
                calcRise, iterations=1):
    """
    General function to calculate the rise time (if calcRise is True) or the set
    time (if calcRise is False) of different kinds of celestial objects.

    :param string objType: Type of celestial object. Must be one of "star",
    sun", "moon", "planet".
    :param float elevation: Horizontal elevation (degree).
    :param float longitude: Geographical longitude of observer (degree).
    :param float latitude: Geographical latitude of observer (degree).
    :param function objPosFunc: This is a function which returns the position
    (ra [h] / dec [degree]) of the celestial object. The function is called with
    a Python datetime object as parameter.
    :param date date: Date (python date object).
    :param float T0: Initial value for rise/set time (hours as float).
    :param boolean calcRise: True to calculate rise time. False to calculate
    set time.
    :param int iterations: Number of iterations.
    :return: Returns a python datetime object.
    """
    baseDate = datetime.datetime.combine(date, datetime.time())
    T = T0

    for i in range(0, iterations):
        localSiderealTime = (sidereal_time(date.year, date.month, date.day, T) +
                             longitude / 15.0)

        objRa, objDec = objPosFunc(baseDate + datetime.timedelta(hours=T))

        hourAngle = localSiderealTime - objRa # hour angle of object (tau)
        if hourAngle > 12:
            hourAngle -= 24
        elif hourAngle < -12:
            hourAngle += 24

        x = ((sin(elevation) - sin(latitude) * sin(objDec)) /
             (cos(latitude) * cos(objDec)))
        if abs(x) > 1:
            return None # no rise/set, always or never visible
        t = acos(x) / 15.0
        if calcRise:
            t = -t

        if objType == "star":
            n = 1.0027379
        elif objType == "sun":
            n = 1
        elif objType == "moon" or objType == "planet":
            if i == 0:
                if objType   == "moon":   n = 1.0027 - 0.0366
                elif objType == "planet": n = 1.0027 - 0
            else:
                n = 1.0027 - ((objRa - objRa_last) / (T - T_last))
        else:
            raise Exception("invalid object type: %s" % objType)

        T_last = T
        objRa_last = objRa
        T = T + ((t - hourAngle) / n)

    delta = datetime.timedelta(hours=T)
    return (baseDate + delta)

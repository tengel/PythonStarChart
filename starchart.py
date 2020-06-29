#!/usr/bin/env python3
"""
Qt star chart GUI.


----
"""

import sys
import os
import math
import datetime
import pathlib
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
import astro
import kepler

class DrawArea(QWidget):
    BORDER=30
    def __init__(self, parent=None):
        super().__init__(parent)
        self.azGridEnabled = True
        self.horizonEnabled = True
        self.size = 600
        self.setMinimumSize(self.size, self.size);
        pal = self.palette()
        pal.setColor(QPalette.Background, Qt.white)
        self.setAutoFillBackground(True)
        self.setPalette(pal)
        self.setSizePolicy(QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding))

    def drawAzGrid(self, painter):
        painter.setPen(Qt.SolidLine)
        painter.setPen(QColor("green"))
        north = QPoint(self.size / 2 + self.BORDER, self.BORDER)
        south = QPoint(self.size / 2 + self.BORDER, self.size + self.BORDER)
        west = QPoint(self.size + self.BORDER, self.size / 2 + self.BORDER)
        east = QPoint(self.BORDER, self.size / 2 + self.BORDER)
        painter.drawLine(north, south);
        painter.drawLine(east, west);
        painter.drawText(north.x() - 5, north.y() - 2, "N")
        painter.drawText(west.x() + 2, west.y() + 5, "W")
        painter.drawText(south.x() - 5, south.y() + 13, "S")
        painter.drawText(east.x() - 13, east.y() + 5, "E")

        r = self.size / 6
        cxy = self.size / 2 + self.BORDER
        painter.drawEllipse(QPoint(cxy, cxy), r, r)
        painter.drawEllipse(QPoint(cxy, cxy), r * 2, r * 2)
        painter.drawEllipse(QPoint(cxy, cxy), r * 3, r * 3)
        painter.drawText(cxy + 2, cxy + r - 2, "60°") # elevation
        painter.drawText(cxy + 2, cxy + 2 * r - 2, "30°")
        painter.drawText(cxy + 2, cxy + 3 * r - 2, "0°")
        painter.drawText(north.x() - 10, north.y() - 15, "180°") # azimuth
        painter.drawText(south.x() - 5, south.y() + 28, "0°")
        painter.rotate(-90)
        painter.drawText(west.y() * -1 - 10, west.x() + 26, "90°")
        painter.drawText(east.y() * -1 - 10, east.x() - 15, "270°")
        painter.rotate(90)

    def horizontal2area(self, azimuth, elevation):
        """
        Convert horizontal coordinates to paint area
        Returns (x, y)
        """
        hpixel = self.size / 2 / 90.0
        xoff = math.sin(math.radians(azimuth)) * (90 - elevation) * hpixel
        yoff = math.cos(math.radians(azimuth)) * (90 - elevation) * hpixel
        x = self.size / 2 + self.BORDER + xoff
        y = self.size / 2 + self.BORDER + yoff
        return (x, y)

    def paintEvent(self, event):
        painter = QPainter(self)
        if self.horizonEnabled:
            painter.setPen(QPen(QBrush(QColor("blue"), Qt.SolidPattern), 2))
            painter.drawEllipse(self.BORDER, self.BORDER, self.size, self.size)
        if self.azGridEnabled:
            self.drawAzGrid(painter)
        for o in self.objects:
            o.draw(self, painter)

    def resizeEvent(self, event):
        self.size = min(self.width(), self.height()) - (2 * self.BORDER)

    def addObjects(self, objects):
        self.objects = objects


class Star:
    def __init__(self):
        self.name = ""
        self.rightAscension = 0
        self.declination = 0
        self.apparentMagnitude = 0
        self.azimuth = 0
        self.elevation = 0

    def __str__(self):
        return ("Star: '%s' HR%4d ra:%fh de:%f° az:%f ele:%f mag:%f" %
                (self.name, self.hr, self.rightAscension, self.declination,
                 self.azimuth, self.elevation, self.apparentMagnitude))

    def draw(self, chart, painter):
        painter.setPen(Qt.SolidLine)
        painter.setPen(QColor("black"))
        painter.setBrush(QColor("black"))
        x, y = chart.horizontal2area(self.azimuth, self.elevation)
        if self.apparentMagnitude < 2:
            painter.drawEllipse(QPoint(x, y), 1.5, 1.5)
        elif self.apparentMagnitude > 4:
            painter.drawPoint(x, y)
        else:
            painter.drawEllipse(QPoint(x, y), 0.9, 0.9)


class EqGrid:
    def __init__(self, engine):
        self.equator = []
        self.decGrid = []
        self.raGrid = []
        self.centerAzimuth, self.centerElevation = engine.equatorial2horizontal(0.0, 90.0)
        self.equator = [engine.equatorial2horizontal(ra, 0) for ra in range(0, 25)]
        for dec in range(-30, 90, 30):
            self.decGrid.append([dec, [engine.equatorial2horizontal(ra, dec)
                                       for ra in range(0, 25)]])
        for ra in range(0, 25, 2):
            self.raGrid.append([ra, [engine.equatorial2horizontal(ra, dec)
                                     for dec in range(-30, 120, 30)]])

    def draw(self, chart, painter):
        painter.setPen(Qt.DashLine)
        painter.setPen(QColor("red"))
        painter.setBrush(Qt.NoBrush)
        x, y = chart.horizontal2area(self.centerAzimuth,
                                    self.centerElevation)
        areaCoords = [chart.horizontal2area(hor[0], hor[1]) for hor in self.equator]
        polygon = QPolygon([QPoint(a[0], a[1]) for a in areaCoords])
        painter.drawPolygon(polygon)
        for g in self.raGrid:
            areaCoords = [chart.horizontal2area(hor[0], hor[1]) for hor in g[1]]
            polygon = QPolygon([QPoint(a[0], a[1]) for a in areaCoords])
            painter.drawPolyline(polygon)
            painter.drawText(areaCoords[0][0] - 10, areaCoords[0][1] + 13, "%d h" % g[0])
        for g in self.decGrid:
            areaCoords = [chart.horizontal2area(hor[0], hor[1]) for hor in g[1]]
            polygon = QPolygon([QPoint(a[0], a[1]) for a in areaCoords])
            painter.drawPolyline(polygon)
            painter.drawText(areaCoords[0][0] - 10, areaCoords[0][1] + 13, "%d°" % g[0])

class Ecliptic:
    def __init__(self, engine):
        self.line = []
        for lon in range(0, 360, 10):
            rightAscension, declination = astro.geoEcl2geoEqua(0.0, lon)
            azimuth, elevation = engine.equatorial2horizontal(rightAscension / 15,
                                                              declination)
            self.line.append((azimuth, elevation))

    def draw(self, chart, painter):
        painter.setPen(QPen(QBrush(QColor("orange"), Qt.SolidPattern), 2))
        painter.setBrush(Qt.NoBrush)
        c = [chart.horizontal2area(x[0], x[1]) for x in self.line]
        polygon = QPolygon([QPoint(x[0], x[1]) for x in c])
        painter.drawPolygon(polygon)


class Planet:
    def __init__(self, engine, p):
        self.p = p
        self.azimuth, self.elevation =  engine.equatorial2horizontal(
            p.ra / 15, p.declination)

    def draw(self, chart, painter):
        painter.setPen(Qt.SolidLine)
        painter.setPen(QColor("red"))
        painter.setBrush(QColor("red"))
        x, y = chart.horizontal2area(self.azimuth, self.elevation)
        painter.drawEllipse(QPoint(x, y), 3, 3)
        painter.drawText(x, y, "  %s" % self.p.name)

    def __str__(self):
        return ("Planet %s ra:%f dec:%f az:%f el:%f" %
              (self.p.name, self.p.ra, self.p.declination,
               self.azimuth, self.elevation))

class Sun:
    def __init__(self, engine):
        raDec = astro.calcPositionSun(astro.julian_date(engine.dateTime))
        self.azimuth, self.elevation = engine.equatorial2horizontal(raDec[0] / 15,
                                                                    raDec[1])
    def draw(self, chart, painter):
        painter.setPen(Qt.SolidLine)
        painter.setPen(QColor("black"))
        painter.setBrush(QColor("Yellow"))
        x, y = chart.horizontal2area(self.azimuth, self.elevation)
        painter.drawEllipse(QPoint(x, y), 6, 6)
        painter.drawText(x, y, "  Sun")

class Moon:
    def __init__(self, engine):
        raDec = astro.calcPositionMoon(astro.julian_date(engine.dateTime))
        self.azimuth, self.elevation = engine.equatorial2horizontal(raDec[0] / 15,
                                                                    raDec[1])
    def draw(self, chart, painter):
        painter.setPen(Qt.SolidLine)
        painter.setPen(QColor("black"))
        painter.setBrush(QColor("LightBlue"))
        x, y = chart.horizontal2area(self.azimuth, self.elevation)
        painter.drawEllipse(QPoint(x, y), 5, 5)
        painter.drawText(x, y, "  Moon")

class Catalog:
    def __init__(self):
        self.list = []
        f = open(pathlib.Path(__file__).parent / "data" / "bs" / "catalog")
        for line in f:
            try:
                HR = line[0:4]
                Name = line[5:14]
                DM = line[14:25]
                HD = line[25:31]
                SAO = line[31:37]
                FK5 = line[37:41]
                IRflag = line[41]
                r_IRflag = line[42]
                Multiple = line[43]
                ADS = line[44:49]
                ADScomp = line[49:51]
                VarID = line[51:60]
                RAh1900 = line[60:62]
                RAm1900 = line[62:64]
                RAs1900 = line[64:68]
                DE_1900 = line[68]
                DEd1900 = line[69:71]
                DEm1900 = line[71:73]
                DEs1900 = line[73:75]
                RAh = int(line[75:77])
                RAm = int(line[77:79])
                RAs = float(line[79:83])
                DE_ = -1 if line[83] == "-" else +1
                DEd = int(line[84:86])
                DEm = int(line[86:88])
                DEs = int(line[88:90])
                GLON = line[90:96]
                GLAT = line[96:102]
                Vmag = float(line[102:107])
                n_Vmag = line[107]
                u_Vmag = line[108]
                B_V = line[109:114]
                u_B_V = line[114]
                U_B = line[115:120]
                u_U_B = line[120]
                R_I = line[121:126]
                n_R_I = line[126]
                SpType = line[127:147]
                n_SpType = line[147]
                pmRA = line[148:154]
                pmDE = line[154:160]
                n_Parallax = line[160]
                Parallax = line[161:166]
                RadVel = line[166:170]
                n_RadVel = line[170:174]
                l_RotVel = line[174:176]
                RotVel = line[176:179]
                u_RotVel = line[179:180]
                Dmag = line[180:184]
                Sep = line[184:190]
                MultID = line[190:194]
                MultCnt = line[194:196]
                NoteFlag = line[196:197]
                s = Star()
                s.name = Name
                s.hr = int(HR)
                s.rightAscension = RAh + (RAm / 60.0) + (RAs / 60.0 / 60.0)
                s.declination = (DEd + (DEm / 60.0) + (DEs / 60.0 / 60.0)) * DE_
                s.apparentMagnitude = Vmag
                self.list.append(s)
            except Exception as e:
                pass
        f.close()


class Engine:
    def __init__(self):
        self.catalog = Catalog()
        self.maxMagnitude = 4
        self.eqGridEnabled = True
        self.eclipticEnabled = True
        self.earth = kepler.Earth()
        self.planets = [kepler.Mercury(), kepler.Venus(), kepler.Mars(),
                        kepler.Jupiter(), kepler.Saturn(), kepler.Uranus(),
                        kepler.Neptune()]

    def setLocation(self, lat, lon):
        self.lat = lat
        self.lon = lon

    def setTime(self, dateTime):
        """
        Set local time with timezone information. Internal calculation is done in UTC.
        """
        self.dateTime = dateTime - dateTime.utcoffset()

    def equatorial2horizontal(self, rightAscension, declination):
        """
        Convert geocentric, equatorial to topocentric, equatorial

        Uses current time and location.

        :param float rightAscension: Geocentric, equatorial right ascension in h
        :param float declination: Geocentric equatorial declination in degree
        :return: (azimut, elevation)
        """
        hourAngle = (self.localSiderealTime - rightAscension) * 15.0
        azimuth, elevation = astro.geoEqua2geoHori(hourAngle, self.lat, declination)
        return (azimuth, elevation)

    def update(self):
        self.objects = []
        utcHour = self.dateTime.hour + self.dateTime.minute / 60.0
        siderealTime = astro.sidereal_time(self.dateTime.year, self.dateTime.month,
                                           self.dateTime.day, utcHour)
        self.localSiderealTime = siderealTime + self.lon / 15.0 # in h
        for star in self.catalog.list:
            if star.apparentMagnitude > self.maxMagnitude:
                continue
            star.azimuth, star.elevation =  self.equatorial2horizontal(
                star.rightAscension, star.declination)
            self.objects.append(star)

        self.earth.calcHeliocentric(self.dateTime)
        for p in self.planets:
            p.calcHeliocentric(self.dateTime)
            p.calcGeocentric(self.earth)
            self.objects.append(Planet(self, p))

        if self.eqGridEnabled:
            self.objects.append(EqGrid(self))
        if self.eclipticEnabled:
            self.objects.append(Ecliptic(self))
        self.objects.append(Sun(self))
        self.objects.append(Moon(self))

    def getObjects(self):
        return self.objects


class Chart(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.chartLayout = QVBoxLayout()
        self.area = DrawArea(self)
        self.chartLayout.addWidget(self.area)
        self.setLayout(self.chartLayout)
        self.engine = Engine()
        self.engine.setLocation(53.14, 8.19)
        self.engine.setTime(datetime.datetime.now().astimezone())
        self.engine.update()
        self.area.addObjects(self.engine.getObjects())

    def update(self):
        self.engine.update()
        self.area.addObjects(self.engine.getObjects())
        self.area.update()

    def showAzGrid(self, t):
        self.area.azGridEnabled = (t == Qt.Checked)
        self.update()

    def showEqGrid(self, t):
        self.engine.eqGridEnabled = (t == Qt.Checked)
        self.update()

    def showHorizon(self, t):
        self.area.horizonEnabled = (t == Qt.Checked)
        self.update()

    def showEcliptic(self, t):
        self.engine.eclipticEnabled = (t == Qt.Checked)
        self.update()

    def showEquator(self):
        pass

    def setMaxMagnitude(self, mag):
        self.engine.maxMagnitude = mag
        self.update()

    def setDateTime(self, t):
        self.engine.setTime(t.toPyDateTime().astimezone())
        self.update()

    def setLocation(self, latitude, longitude):
        self.engine.setLocation(latitude, longitude)
        self.update()

if __name__ == "__main__":
    app = QApplication([])

    window = QWidget()
    window.setWindowTitle(os.path.basename(sys.argv[0]))

    quitShortcut = QShortcut(QKeySequence(Qt.Key_Q), window);
    quitShortcut.activated.connect(window.close)

    chart = Chart()
    timeLabel = QLabel("Local Time:")
    dateTime = QDateTimeEdit(QDateTime.currentDateTime())
    dateTime.setDisplayFormat("yyyy-MM-dd\t HH:mm")
    dateTime.dateTimeChanged.connect(chart.setDateTime)
    locBox = QGroupBox("Location:")
    latBox = QDoubleSpinBox()
    latBox.setMinimum(-90)
    latBox.setMaximum(90)
    latBox.setValue(53.14)
    latBox.valueChanged.connect(lambda: chart.setLocation(latBox.value(),
                                                          lonBox.value()))
    lonBox = QDoubleSpinBox()
    lonBox.setValue(8.19)
    lonBox.setMinimum(-180)
    lonBox.setMaximum(180)
    lonBox.valueChanged.connect(lambda: chart.setLocation(latBox.value(),
                                                          lonBox.value()))
    magLabel = QLabel("Max. apparent magnitude:")
    magBox = QSpinBox()
    magBox.setMaximum(8)
    magBox.setValue(chart.engine.maxMagnitude)
    magBox.valueChanged.connect(chart.setMaxMagnitude)
    optionsLabel = QLabel("Display options:")
    horizonCb = QCheckBox("Horizon")
    horizonCb.setCheckState(Qt.Checked if chart.area.horizonEnabled else
                            Qt.Unchecked)
    horizonCb.stateChanged.connect(chart.showHorizon)
    equatorCb = QCheckBox("Equator")
    eclipticCb = QCheckBox("Ecliptic")
    eclipticCb.setCheckState(Qt.Checked if chart.engine.eclipticEnabled else
                             Qt.Unchecked)
    eclipticCb.stateChanged.connect(chart.showEcliptic)
    eqGridCb = QCheckBox("Equatorial grid")
    eqGridCb.setCheckState(Qt.Checked if chart.engine.eqGridEnabled else
                           Qt.Unchecked)
    eqGridCb.stateChanged.connect(chart.showEqGrid)
    azGridCb = QCheckBox("Azimuthal grid")
    azGridCb.setCheckState(Qt.Checked if chart.area.azGridEnabled else
                           Qt.Unchecked)
    azGridCb.stateChanged.connect(chart.showAzGrid)
    exitBtn = QPushButton("Exit")
    exitBtn.clicked.connect(window.close)

    locLayout = QGridLayout()
    locLayout.addWidget(QLabel("Latitude:"), 0, 0)
    locLayout.addWidget(latBox, 0, 1)
    locLayout.addWidget(QLabel("Longitude:"), 1, 0)
    locLayout.addWidget(lonBox, 1, 1)
    locBox.setLayout(locLayout)

    leftLayout = QVBoxLayout()
    leftLayout.addWidget(timeLabel)
    leftLayout.addWidget(dateTime)
    leftLayout.addWidget(locBox)
    leftLayout.addWidget(magLabel)
    leftLayout.addWidget(magBox)
    leftLayout.addStretch()
    leftLayout.addWidget(optionsLabel)
    leftLayout.addWidget(horizonCb)
    leftLayout.addWidget(equatorCb)
    leftLayout.addWidget(eclipticCb)
    leftLayout.addWidget(eqGridCb)
    leftLayout.addWidget(azGridCb)

    leftLayout.addWidget(exitBtn)

    windowLayout = QHBoxLayout()
    windowLayout.addLayout(leftLayout)
    windowLayout.addWidget(chart)

    window.setLayout(windowLayout)
    window.show()

    app.exec_()

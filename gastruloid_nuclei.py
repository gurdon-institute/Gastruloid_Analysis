# Segment nuclei in DAPI channel excluding aggregates and measure intensity in each channel.
# by Richard Butler, Gurdon Institute Imaging Facility, University of Cambridge
#
# Copyright 2024 Richard Butler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>

import math as maths

from ij import IJ, WindowManager, ImagePlus, ImageStack
from ij.gui import ShapeRoi, TextRoi, Overlay
from ij.process import Blitter, ImageProcessor, ShortProcessor, AutoThresholder, FloodFiller
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.measure import ResultsTable

from java.awt import Color, Font, Canvas


LABELFONT = Font(Font.SANS_SERIF, Font.PLAIN, 12)
fm = Canvas().getFontMetrics(LABELFONT)

sigma = 1.0		# µm
k = 5.0

tolerance = 0.9

minA = 120.0 # µm²
maxA = 600.0


def watershed(ip):
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, tolerance, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)


def fillHoles(ip):
	width = ip.getWidth()
	height = ip.getHeight()
	ff = FloodFiller(ip)
	ip.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if ip.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if ip.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if ip.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if ip.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if ip.get(i)==127:
		    ip.set(i, 0)
		else:
		    ip.set(i, 255)


def onEdge(roi, imp):
	width = imp.getWidth()
	height = imp.getHeight()
	for xy in range(max(width,height)):
		if roi.contains(0,xy) or roi.contains(xy,0) or roi.contains(width-1,xy) or roi.contains(xy,height-1):
			return True
	return False


imp = IJ.getImage()
cal = imp.getCalibration()
minApx = minA / (cal.pixelWidth * cal.pixelHeight)
maxApx = maxA / (cal.pixelWidth * cal.pixelHeight)
stack = imp.getStack()
ol = Overlay()

nucmap = stack.getProcessor(1).duplicate()
sub = nucmap.duplicate()
nucmap.blurGaussian(sigma/cal.pixelWidth)
sub.blurGaussian(k*sigma/cal.pixelWidth)
nucmap.copyBits(sub, 0,0, Blitter.SUBTRACT)
hist = nucmap.getStatistics().getHistogram()
hist = [int(n) for n in hist]
thresh = AutoThresholder().getThreshold( AutoThresholder.Method.Otsu, hist )
nucmap.threshold(int(thresh))
mask = nucmap.convertToByte(False)
fillHoles(mask)
watershed(mask)
mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
composite = ThresholdToSelection().convert(mask)
rois = ShapeRoi(composite).getRois()
nuclei = []
for roi in rois:
	if minApx <= roi.getStatistics().area <= maxApx and not onEdge(roi, imp):
		nuclei.append(roi)

rt = ResultsTable.getResultsTable()
rt.showRowNumbers(False)
for i,nuc in enumerate(nuclei):
	nuc.setStrokeColor(Color.MAGENTA)
	ol.add(nuc)
	row = rt.getCounter()
	rt.setValue("Image", row, imp.getTitle())
	statsA = nuc.getStatistics()
	nucX = statsA.xCentroid * cal.pixelWidth
	nucY = statsA.yCentroid * cal.pixelHeight

	rt.setValue("Nucleus", row, i)
	
	rect = nuc.getBounds()
	txt = str(i)
	labelX = int(rect.x + rect.width/2. - fm.stringWidth(txt)/2.)
	labelY = rect.y - fm.getHeight()
	label = TextRoi(labelX, labelY, txt, LABELFONT)
	label.setStrokeColor(Color.CYAN)
	ol.add(label)
	
	rt.setValue("X", row, nucX)
	rt.setValue("Y", row, nucY)
	rt.setValue("Area", row, statsA.area * cal.pixelWidth * cal.pixelHeight)
	for c in range(2,imp.getNChannels()+1):
		cip = stack.getProcessor(c)
		cip.setRoi(nuc)
		stats = cip.getStatistics()
		rt.setValue("C"+str(c)+" Mean", row, stats.mean)
		rt.setValue("C"+str(c)+" StdDev", row, stats.stdDev)

imp.setOverlay(ol)

rt.show("Results")

# Segment an organoid in the DAPI channel, measure proportion of volume positive for signal in C2, C3, C4 and C2+C3
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
import re
from collections import OrderedDict

from ij import IJ, ImagePlus, ImageStack
from ij.gui import ShapeRoi, Overlay
from ij.process import Blitter, ImageProcessor, ShortProcessor, AutoThresholder, StackStatistics
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.plugin import GaussianBlur3D
from ij.measure import ResultsTable

from java.awt import Color, Font, Canvas


COLOURS = {
    "C1": [Color.RED, "Red"],
    "C2": [Color.GREEN, "Green"],
    "C3": [Color.CYAN, "Cyan"],
    "C4": [Color(255,128,0), "Orange"],
    "C2+C3": [Color(128,0,128), "Purple"]
}
SIGMA = 3.0 # µm
METHODS = {
    1: AutoThresholder.Method.Triangle,
    2: AutoThresholder.Method.MaxEntropy,
    3: AutoThresholder.Method.Otsu,
    4: AutoThresholder.Method.Triangle
}
CELLR = 5.0	# half distance between nucleus centroids, include non-nuclear volume

def getMask(imp, channel, sigma, method):
	cal = imp.getCalibration()
	
	sigmaXY = sigma/cal.pixelWidth
	sigmaZ = sigma/cal.pixelDepth
	
	procstack = ImageStack()
	for z in range(1,nSlices+1):
		procstack.addSlice(stack.getProcessor(imp.getStackIndex(channel,z,1)).duplicate())
	procimp = ImagePlus("proc_c"+str(channel),procstack)
	GaussianBlur3D.blur(procimp, sigmaXY, sigmaXY, sigmaZ)
	
	if channel==2:
		procstack = procimp.getStack()
		for z in range(1,imp.getNSlices()+1):
			procstack.getProcessor(z).log()
	
	stats = StackStatistics(procimp)
	hist = [int(h) for h in stats.getHistogram()]
	#hist[0] = 0
	thresh = AutoThresholder().getThreshold( method, hist )
	thresh = (thresh/float(255) * (stats.max-stats.min) + stats.min)
	out = ImageStack()
	for z in range(1,imp.getNSlices()+1):
		ip = procstack.getProcessor(z)
		ip.threshold(int(thresh))
		mask = ip.convertToByte(False)
		out.addSlice(mask)
	#procimp.flush()
	return out


def getRois(maskstack):
	rois = []
	tts = ThresholdToSelection()
	for z in range(1,maskstack.size()+1):
		mip = maskstack.getProcessor(z)
		mip.setThreshold(255,255, ImageProcessor.NO_LUT_UPDATE)
		composite = tts.convert(mip)
		if composite is None:
			continue
		composite.setPosition(-1,z,1)
		rois.append(composite)
	return rois


imp = IJ.getImage()
nChannels = imp.getNChannels()
nSlices = imp.getNSlices()
cal = imp.getCalibration()
stack = imp.getStack()
ol = Overlay()
rt = ResultsTable.getResultsTable()
rt.showRowNumbers(False)



volmask = OrderedDict()
for c in range(1,nChannels+1):
	volmask["C"+str(c)] = getMask(imp, c, SIGMA, METHODS[c])

doublemask = ImageStack()
for z in range(1,nSlices+1):
	c2mask = volmask["C2"].getProcessor(z)
	c3mask = volmask["C3"].getProcessor(z)
	dip = c2mask.duplicate()
	dip.copyBits(c3mask, 0,0, Blitter.AND)
	doublemask.addSlice(dip)
volmask["C2+C3"] = doublemask

title = imp.getTitle()
tre = re.compile("_(\d+)_DAPI")
m = tre.search(title)
titlen = m.group(1) if m is not None else "?"
cellV = (4/3.)*maths.pi*CELLR**3
volume = {}
for ckey in volmask.keys():
	rois = getRois(volmask[ckey])
	colour,cname = COLOURS[ckey]
	volume[ckey] = 0
	for roi in rois:
		volume[ckey] += roi.getStatistics().area * cal.pixelWidth * cal.pixelHeight * cal.pixelDepth
		roi.setStrokeColor(colour)
		ol.add(roi)
	row = rt.getCounter()
	rt.setValue("Image",row,"hG-"+titlen)
	rt.setValue("Marker",row,ckey)
	rt.setValue("Colour",row,cname)
	rt.setValue("+ve Volume ("+cal.getUnit()+u"\u00b3"+")",row,volume[ckey])
	rt.setValue("+ve Proportion",row,volume[ckey]/volume["C1"])
	rt.setValue("Estimated Cell Count",row,volume[ckey]/cellV)

imp.setOverlay(ol)

rt.show("Results")

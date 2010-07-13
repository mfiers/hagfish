import os
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
import matplotlib as mpl
import matplotlib.transforms
import logging
import optparse

################################################################################
## GLOBALS
##
COLOR1 = '#28A62E'
COLOR2 = '#DF1818'
COLOR3 = '#3D358B'
COLOR4 = '#F0D518'
COLOR5 = '#01BEAF'

################################################################################
## set up logging
##
def getLogger(what, level):
    l = logging.getLogger(what)
    handler = logging.StreamHandler()
    logmark = chr(27) + '[0;37;44m' + what + \
              chr(27) + '[0m ' 
    formatter = logging.Formatter(
        logmark + '%(levelname)-6s %(message)s')
    handler.setFormatter(formatter)
    l.addHandler(handler)
    
    if level >= 2:
        l.setLevel(logging.DEBUG)
    elif level == 1:
        l.setLevel(logging.INFO)
    else:
        l.setLevel(logging.WARNING)

    return l

################################################################################
## Define optparse parameters
##
def getHagfishOptparser():
    parser = optparse.OptionParser()
    parser.add_option('-v', dest='verbose', action="count", 
                      help='Show debug information')
    return parser

def addPlotParameters(parser):
    parser.set_defaults(ntPerBand=1e6)
    parser.add_option('-n', dest='ntPerBand',
                      help='no nucleotides per band')
    parser.add_option('-i', dest='inputFile',
                      help='input file with the coverage data (npz, if not specified, '+
                      'the input file name will be inferred from the sequence Id')
    parser.add_option('-l', dest='log', action='store_true',
                      help='take the log10 of the coverage')

    parser.set_defaults(imageWidth=1000)
    parser.set_defaults(bandHeight=200)
    parser.add_option('-W', dest='imageWidth', help='imageWidth (in px)')
    parser.add_option('-H', dest='bandHeight', help='bandHeight (in px)')

    parser.set_defaults(yfrac=0.98)
    parser.add_option('-Y', dest='yfrac', type='float', help='percentage of the plotted'
                      'fraction that must fall inside the Y boundaries of the graph - use'
                      'this to scale the y axis')
    parser.add_option('--ymax', dest='ymax',
                      help='Alternatively, set a max value for the y axis')

    parser.add_option('-s', dest='start',
                      help='Start position (nt) of the plot')
    parser.add_option('-e', dest='stop',
                      help='Stop position (nt) of the plot')


    parser.add_option('-o', dest='outfile',
                      help='Output file name')

    parser.set_defaults(format=['png'])
    parser.add_option('-f', dest='format', action='append',
                      help='Output format (may be defined more than once)')


class hagfishData:
    
    def __init__(self, options, args):

        self.seqId = args[0]
        self.options = options

        self.l = getLogger('data', options.verbose)
        self.l.info("Loading sequence: %s" % self.seqId)

        if options.inputFile:
            self.inputFile = options.inputFile
        else:
            self.inputFile = os.path.join(
                'combined', 
                '%s.combined.coverage.npz' % self.seqId)

        self.l.info('loading %s' % self.inputFile)
        self.vectors = ['ok', 'high', 'low']
        self.data = np.load(self.inputFile)
        self.ok = self.data['r_ok']
        self.high = self.data['r_high']
        self.low = self.data['r_low']
        #self.low_binned = self.data['r_low_binned']
        #self.high_binned = self.data['r_high_binned']
        self.all = self.low + self.ok + self.high

        self.seqLen = len(self.data['r_ok'])
        self.l.debug('loaded coverage plots of len %d' % self.seqLen) 

        #some stats
        self.okh = self.ok / 2
        self.vectors.append('okh')
        self.medianH = np.median(self.okh)
        self.median = np.median(self.ok)

        #some stuff for plotting
        self.x = np.arange(0, self.seqLen, dtype="int")
        #self.z = np.zeros_like(self.x)
        self.vectors.append('x')

class hagfishPlot:

    def __init__(self, options, data):
        self.l = getLogger('plot', options.verbose)
        self.options = options
        self.data = data

        self.start = 0
        if options.start:
            self.start = int(float(options.start))

        self.stop = self.data.seqLen
        if options.stop:
            self.stop = int(float(options.stop))
            if self.stop > data.seqLen:
                self.stop = data.seqLen
                
        self.plotLen = self.stop - self.start
        self.ntPerBand = int(float(options.ntPerBand))
        if self.ntPerBand > self.plotLen:
            self.ntPerBand = self.plotLen

        self.l.info("Ploting %d nucleotides per band" % self.ntPerBand)
        self.l.info("Plotting from %d to %d (%d nt)" % (
            self.start, self.stop, self.plotLen))

        self.noBands = int(math.ceil(self.plotLen / float(self.ntPerBand)))
        self.l.info("Plotting %d bands" % self.noBands)

        self.yTicks = []
        self.yTicks2 = []
        self.yTickLabels = []
        self.yTickLabels2 = []

        self.tminY, self.tmaxY = 0, 0

        if options.ymax:
            self.maxY = int(options.ymax)
            self.minY = -maxY
        else:
            self.maxY = int(quant(self.data.okh, options.yfrac))
            self.minY = -self.maxY

        self.YCorrPerBand = 2 * self.maxY

        """
        Prepare matplotlib
        """
        mpl.rcParams['axes.labelsize'] = 'x-small'
        mpl.rcParams['xtick.labelsize'] = 'xx-small'
        mpl.rcParams['ytick.labelsize'] = 'xx-small'

        self.figWidth = self.options.imageWidth / 100.0
        self.bandHeight = self.options.bandHeight / 100.0
        self.figHeight = self.bandHeight * self.noBands
        self.fig = plt.figure(figsize=(self.figWidth, self.figHeight))
        figtb = 0.2 / self.noBands
        self.l.debug("no bands: %d" % self.noBands)
        self.axCoords = 0.05, figtb, 0.9, 1 - (2*figtb)
        self.l.debug("creating axes with coords lbrt %.4f %.4f %.4f %.4f " %
                     self.axCoords)

        self.ax = self.fig.add_axes(self.axCoords)

        if self.start == 0:
            self.ax.set_title('Coverage plot for "%s" (0 to %.0e)' % (
                self.data.seqId, self.stop))
        else:
            self.ax.set_title('Coverage plot for "%s" (%.0e to %.0e)' % (
                self.data.seqId, self.start, self.stop))

    def plotBands(self, bandPlotter, layers=1):
        self.bands = []
        yCorrect = 0

        for band in range(self.noBands):
            b = bandPlotter(self, band, self.options)
            yCorrect -= self.YCorrPerBand
            b.setYcorrection(yCorrect)

            self.yTicks.append(yCorrect)
            self.yTickLabels.append("%.0e" % b.start)

            b.setYticks2()
            self.bands.append(b)
            if band == 0:
                self.tmaxY = yCorrect + self.maxY
            if band == self.noBands -1:
                self.tminY = yCorrect - self.maxY

            b.plotBand()
            
        y = self.ax.get_yaxis()
        y.set_ticks(self.yTicks)
        y.set_ticklabels(self.yTickLabels)
        self.ax.set_ylim(self.tminY,self.tmaxY)

        y2 = self.ax.twinx()
        y2.get_yaxis().set_ticks(self.yTicks2)
        y2.get_yaxis().set_ticklabels(self.yTickLabels2)
        y2.set_ylim(self.tminY,self.tmaxY)

    def save(self):
        self.ax.set_xlim(0, self.ntPerBand)
        if self.options.outfile:            
            outFileName = self.options.outfile
        else:
            outFileName = self.data.seqId

        if self.options.start or self.options.stop:
            outFileName += "_%d_%d" % (self.start, self.stop)

        for f in self.options.format:
            plt.savefig('%s.%s' % (outFileName, f), dpi=100)


class hagfishPlotBand:
    
    def __init__(self, plot, band, options):

        self.l = getLogger('band.%03d' % band, options.verbose)
        self.plot = plot
        self.ax = plot.ax
        self.data = self.plot.data
        self.band = band
        self.start = self.plot.start + (band * self.plot.ntPerBand)
        self.stop = self.start + (self.plot.ntPerBand - 1)
        self.l.debug(
            'printing from %d to %d' % (
                self.start, self.stop))

        for v in self.data.vectors:
            self.__dict__[v] = self.data.__dict__[v][self.start:self.stop]
            
        #special form of x
        self.locx = self.x - self.start
        
    def plot(self, layer):
        """
        Should be overridden
        """
        pass
    
    def setYticks(self):
        self.plot.yTicks2.append(thisBandCorr - hmed)
        self.plot.yTicks2.append(thisBandCorr)
        self.plot.yTicks2.append(thisBandCorr + hmed)

        self.plot.yTickLabels2.append("%d" % -fhmed)
        self.plot.yTickLabels2.append("%d" % 0)
        self.plot.yTickLabels2.append("%d" % hmed)

    

    def setYcorrection(self, value):
        self.yCorrection = value
        self.bandTop = self.yCorrection + \
                       (0.5 * self.plot.YCorrPerBand)
        self.bandBottom = self.yCorrection - \
                       (0.5 * self.plot.YCorrPerBand)
        self.l.debug("defining band stats")
        self.l.debug("   ycorrection %s" % self.yCorrection)
        self.l.debug("   bandTop     %s" % self.bandTop)
        self.l.debug("   bandBottom  %s" % self.bandBottom)
        #TransformedBbox(rect[0].get_bbox(), ax.transData) 
        #self.bbox = matplotlib.transforms.TransformedBbox(
        #     matplotlib.transforms.Bbox(
        #        np.array([[0,self.bandBottom],
        #                  [self.plot.ntPerBand,self.bandTop]])),
        #    self.ax.transData)
        _pd = [
            (mpath.Path.MOVETO, (0, self.bandTop)),
            (mpath.Path.LINETO, (self.plot.ntPerBand, self.bandTop)),
            (mpath.Path.LINETO, (self.plot.ntPerBand, self.bandBottom)),
            (mpath.Path.LINETO, (0, self.bandBottom)),
            (mpath.Path.CLOSEPOLY, (0, self.bandTop)),
            ]
        codes, verts = zip(*_pd)
        self.bpath = mpath.Path(verts, codes)
        self.transform = self.ax.get_transform()
        self.patch = mpatches.PathPatch(self.bpath, facecolor='orange', alpha=0, lw=0)
        self.ax.add_patch(self.patch)



################################################################################
## extra routines for dealing with numpy arrays
##
    
def quant(x, w):
    xs = np.sort(x)
    return xs[int(len(xs) * w)]

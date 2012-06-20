import os
import sys
import math
import copy

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
import matplotlib as mpl
import matplotlib.transforms
import logging
import optparse

from hagfish_file_util import *

################################################################################
## GLOBALS
##
COLOR1 = '#2C722F' #green
COLOR2 = '#8C201E' #red 
COLOR3 = '#224E73' #purple
COLOR4 = '#DEA000' #yellow
COLOR5 = '#078C7A' #cyan
COLOR6 = '#EE1904' #lighter red
COLOR7 = '#32AE1C' # lighter green

COLDARKGREEN = '#007D1C'
COLDARKRED = '#A61000'
COLDARKBLUE = '#015666'
COLDARKYELLOW = '#A65100'

COLLIGHTGREEN = '#38E05D'
COLLIGHTRED = '#FF5240'
COLLIGHTBLUE = '#37B6CE'
COLLIGHTYELLOW = '#FF9D40'

COLDARKGREEN = '#007D1C'
COLDARKRED = '#9F0025'
COLDARKPURPLE = '#42036F'
COLDARKBLUE = '#06266F'
COLDARKYELLOW = '#A66F00'

COLLIGHTGREEN = '#38E05D'
COLLIGHTRED = '#B12D4C'
COLLIGHTBLUE = '#4671D5'
COLLIGHTYELLOW = '#FFBF40'
COLLIGHTPURPLE = '#963FD5'

COLLIST = [COLDARKGREEN, COLDARKRED, COLDARKPURPLE, COLDARKBLUE, COLDARKYELLOW]

COLMAP1 = mpl.colors.LinearSegmentedColormap.from_list(
    'COLMAP1', 
    [COLLIGHTPURPLE, COLLIGHTRED , COLLIGHTRED ,
     'black', COLDARKGREEN, COLDARKGREEN,
     COLLIGHTGREEN ], N=200)

COLMAP2 = mpl.colors.LinearSegmentedColormap.from_list(
    'COLMAP1', 
    [COLLIGHTRED , COLDARKRED, 'black', COLDARKGREEN, COLLIGHTGREEN ], N=200)


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


def addBasePlotParameters(parser):
    parser.set_defaults(ntPerBand=-1)

    parser.add_option('-n', dest='ntPerBand',
                      help='no nucleotides per band')
    
    parser.set_defaults(imageWidth=1000)
    parser.set_defaults(bandHeight=200)
    parser.add_option('-W', dest='imageWidth', type='int', help='imageWidth (in px)')
    parser.add_option('-H', dest='bandHeight', type='int', help='bandHeight (in px)')

    parser.add_option('-s', dest='start',
                      help='Start position (nt) of the plot')
    parser.add_option('-e', dest='stop',
                      help='Stop position (nt) of the plot')

    parser.add_option('-l', dest='library', action='append',
                      help='Library to load - omit to load all')

    parser.add_option('-o', dest='outfile',
                      help='Output file name')

    parser.set_defaults(format=['png'])
    parser.add_option('-f', dest='format', action='append',
                      help='Output format (may be defined more than once)')

    parser.set_defaults(format=['png'], dpi=100)
    parser.add_option('--dpi', dest='dpi', type='int',
                      help='dpi of the image, pixel calculations are based on dpi 100, setting dpi to 200 will double the x/y pixel size of your image)')


def addPlotParameters(parser):
    addBasePlotParameters(parser)
    parser.add_option('-i', dest='inputFile',
                      help='input file with the coverage data (npz, if not specified, '+
                      'the input file name will be inferred from the sequence Id')

    parser.set_defaults(yfrac=0.98)
    parser.add_option('-Y', dest='yfrac', type='float', help='percentage of the plotted'
                      'fraction that must fall inside the Y boundaries of the graph - use'
                      'this to scale the y axis')
    parser.add_option('--ymax', dest='ymax',
                      help='Alternatively, set a max value for the y axis')
    parser.add_option('-Q', dest='quick', action='store_true',
                      help='Create a "light" version of this graph (if implemented)')

def addBinPlotParameters(parser):
    addBasePlotParameters(parser)


def loadLib(o, base, vectors):
    l.info('loading library from %s' % base)
    
    for vc in vectors:
        d = np_load(base, 'r_' + vc)
        if not o.__dict__.has_key(vc):
            o.__dict__[vc] = d
            o.data[vc] = o.__dict__[vc]
        else:
            o.__dict__[vc] += d


    # self.ok = np_load(file_base, 'r_ok')
    # self.high = np_load(file_base, 'r_high')
    # self.low = np_load(file_base, 'r_low')
    
    # self.ok_ends = np_load(file_base, 'r_ok_ends')
    # self.high_ends = np_load(file_base, 'r_high_ends')
    # self.low_ends = np_load(file_base, 'r_low_ends')
    
    # self.data['ok'] = self.ok
    # self.data['low'] = self.low
    # self.data['high'] = self.high
    
    # self.data['ok_ends'] = self.ok_ends
    # self.data['low_ends'] = self.low_ends
    # self.data['high_ends'] = self.high_ends

class hagfishData:
    
    def __init__(self, options, args, seqId = None, inputFile=None, libraries=None):

        if seqId:
            self.seqId = seqId
        else:
            self.seqId = args[0]

        self.options = options

        self.l = getLogger('data', options.verbose)
        self.l.info("Loading sequence: %s" % self.seqId)

        self.vectors = ['ok', 'high', 'low', 
                        'ok_ends', 'high_ends', 'low_ends']

        #self.data = np.load(self.inputFile)
        self.data = {}

        lib2load = []
        if libraries:
            lib2load = libraries
        elif options.library:
            lib2load = options.library

        if lib2load:
            for i, lib in enumerate(lib2load):
                self.l.info("Loading library %s" % lib)
                file_base= os.path.join(
                    'coverage', lib, 
                    '%s.coverage' % self.seqId)
                loadLib(self, file_base, self.vectors)
        else:
            self.l.info("Loading combined library")
            file_base = os.path.join(
                'combined', 
                '%s' % self.seqId)
            loadLib(self, file_base, self.vectors)

        #see if there is gapdata to load
        gapbase = os.path.join('gaps', self.seqId)
        if np_exists(gapbase, 'nns'):
            self.l.info("loading gap data")
            self.nns = np_load(gapbase, 'nns').astype(float)
            self.vectors.append('nns')
        else:
            self.l.info("Cannot find gap data")
            self.nns = None
        self.data['nns'] = self.nns
            
        self.all = self.low + self.ok + self.high

        self.seqLen = len(self.ok)
        self.l.debug('loaded coverage plots of len %d' % self.seqLen) 

        #some stats
        self.okh = self.ok / 2
        self.okh_ends = self.ok_ends / 2
        self.vectors.append('okh')
        self.vectors.append('okh_ends')
        self.medianH = np.median(self.okh)
        self.median = np.median(self.ok)
        self.average = np.average(self.ok)
        self.min = np.min(self.ok)
        self.max = np.max(self.ok)
        self.l.info("stats: median ok is %s (h %s)" % (self.median, self.medianH))
        self.l.info("stats: average ok is %s" % (self.average))
        self.l.info("stats: min, max ok is %s, %s" % (self.min, self.max))
        #some stuff for plotting
        self.x = np.arange(0, self.seqLen, dtype="int")
        #self.z = np.zeros_like(self.x)
        self.vectors.append('x')


class hagfishPlot:

    def __init__(self, options, data, title=None, data2=None, ymax=None, tag=""):
        
        self.tag = tag
        self.l = getLogger('plot', options.verbose)
        self.options = options

        self.data = data
        self.data2 = data2

        self.start = 0
        if options.start:
            self.start = int(float(options.start))
        if self.start > data.seqLen:
            self.l.critical("start plot is greater than sequence length")
            sys.exit(-1)

        self.stop = self.data.seqLen
        if options.stop:
            self.stop = int(float(options.stop))
            if self.stop > data.seqLen:
                self.stop = data.seqLen
            self.l.info("plot stops at %d" % self.stop)

        if self.start > self.stop:
            self.l.critical("plot start is greater than plot stop")
            sys.exit(-1)
                
        self.plotLen = self.stop - self.start
        self.ntPerBand = int(float(options.ntPerBand))

        if self.ntPerBand == -1:
            self.l.debug("nt per band is not specified")                    
            self.ntPerBand = int(1e6)            
            self.noBands = int(math.ceil(self.plotLen / float(self.ntPerBand)))
            while self.noBands > 5:
                self.ntPerBand += 1e6
                self.noBands = int(math.ceil(self.plotLen / float(self.ntPerBand)))
        else:
            self.noBands = int(math.ceil(self.plotLen / float(self.ntPerBand)))

            
        if self.ntPerBand > self.plotLen:
            self.ntPerBand = self.plotLen

        self.l.info("Ploting %d nucleotides per band" % self.ntPerBand)
        self.l.info("Plotting from %d to %d (%d nt)" % (
            self.start, self.stop, self.plotLen))


        self.l.info("Plotting %d bands" % self.noBands)

        self.yTicks = []
        self.yTicks2 = []
        self.yTickLabels = []
        self.yTickLabels2 = []

        self.tminY, self.tmaxY = 0, 0

        if ymax:
            self.l.info("setting ymax (func def) to %s" % ymax)
            self.maxY = ymax
            self.minY = -self.maxY
        elif options.ymax:
            self.l.info("setting ymax (--ymax) to %s" % options.ymax)
            self.maxY = int(options.ymax) / 2
            self.minY = -self.maxY
        else:
            self.maxY = int(quant(self.data.okh, options.yfrac))
            self.minY = -self.maxY
            self.l.info("setting ymax (calc) to %s" % self.maxY)

        self.YCorrPerBand = 2 * self.maxY

        """
        Prepare matplotlib
        """
        mpl.rcParams['axes.labelsize'] = 'x-small'
        mpl.rcParams['axes.titlesize'] = 'small'
        mpl.rcParams['xtick.labelsize'] = 'xx-small'
        mpl.rcParams['ytick.labelsize'] = 'xx-small'

        self.figWidth = self.options.imageWidth / 100.0
        self.bandHeight = float(self.options.bandHeight) / 100.0
        tbh = self.bandHeight * self.noBands
        self.figHeight = tbh + 1
        self.fig = plt.figure(figsize=(self.figWidth, self.figHeight))
        figtb = 1.0 / (3 * (tbh+1))
        self.l.debug("no bands: %d " % self.noBands)
        self.l.debug("top/bottom margin: %f (0.5 / %f ) " % (figtb, tbh))
        #lbrt
        self.axCoords = 0.08, figtb, 0.85, 1 - (2*figtb)
        self.l.debug("creating axes with coords lbrt %.4f %.4f %.4f %.4f " %
                     self.axCoords)

        self.ax = self.fig.add_axes(self.axCoords)

        if not title:
            if self.tag:
                title = '%s plot "%s"' % (self.tag, self.data.seqId)
            else:
                title = '"%s"' % self.data.seqId
            if options.library:
                title += " (%s)" % ", ".join(options.library)
            if self.start == 0:
                title += " (0 to %.2e)" % (self.stop)
            else:
                title += " (%.2e to %.2e)" % (self.start, self.stop)

            self.ax.set_title(title)
        else:
            self.ax.set_title(title)

    def plotBands(self, bandPlotter):
        yCorrect = 0

        for band in range(self.noBands):
            b = bandPlotter(self, band, self.options)
            yCorrect -= self.YCorrPerBand
            b.setYcorrection(yCorrect)
            self.l.debug("start plotting band %d (corr %s)" % (band, yCorrect))

            self.yTicks.append(yCorrect)
            self.yTickLabels.append("%.2e" % b.start)

            b.setYticks2()
            if band == 0:
                self.tmaxY = yCorrect + self.maxY
            if band == self.noBands -1:
                self.tminY = yCorrect - self.maxY
            b.plotBand()
            del b
                
            
        y = self.ax.get_yaxis()
        y.set_ticks(self.yTicks)
        y.set_ticklabels(self.yTickLabels)
        self.ax.set_ylim(self.tminY,self.tmaxY)

        y2 = self.ax.twinx()
        y2.get_yaxis().set_ticks(self.yTicks2)
        y2.get_yaxis().set_ticklabels(self.yTickLabels2)
        y2.set_ylim(self.tminY,self.tmaxY)

    def save(self, libname="", tag=""):

        if tag:
            usetag = tag
        elif self.tag:
            usetag = self.tag
        else:
            usetag = ""

        self.ax.set_xlim(0, self.ntPerBand)

        if self.options.outfile:            
            outFileName = self.options.outfile
        else:
            outFileName = self.data.seqId

        if self.options.start or self.options.stop:
            outFileName += "_%d_%d" % (self.start, self.stop)

        if libname:
            outFileName += "_%s" % libname
        elif self.options.library:
                outFileName += '_%s' % "".join(self.options.library)

        if usetag:
            outFileName += '_%s' % usetag
            
        #plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

        for f in self.options.format:
            self.l.info("writing to %s.%s" % (outFileName, f))
            plt.savefig('%s.%s' % (outFileName, f), dpi=self.options.dpi, 
                        bbox_inches='tight')


class Dummy:
    pass

class hagfishPlotBand:
    
    def __init__(self, plot, band, options):

        self.l = getLogger('band.%03d' % band, options.verbose)
        self.plot = plot
        self.options = options
        self.ax = plot.ax
        self.data = self.plot.data
        self.data2 = self.plot.data2

        self.band = band
        self.start = self.plot.start + (band * self.plot.ntPerBand)
        self.stop = self.start + (self.plot.ntPerBand - 1)
        self.l.debug( 'printing from %d to %d' % ( self.start, self.stop))

        self.l.debug("preprocessing band data")

        self.d1 = Dummy()
        if self.data2: self.d2 = Dummy()

        for v in self.data.vectors:
            setattr(self, v, self.data.__dict__[v][self.start:self.stop])
            setattr(self.d1, v, getattr(self, v))
            if self.data2:
                setattr(self.d2, v, self.data2.__dict__[v][self.start:self.stop])
                    
        #special form of x
        self.locx = self.x - self.start
        self.zero = np.zeros_like(self.locx)

    def plot(self):
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
        self.l.debug("defining band from %s to %s" % (self.bandTop, self.bandBottom))
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

import gzip
import cPickle
import logging

l = logging.getLogger('hagfish')
handler = logging.StreamHandler()
logmark = chr(27) + '[0;37;44mHAGFISH' + \
          chr(27) + '[0m ' 

formatter = logging.Formatter(
    logmark + '%(levelname)-6s %(message)s')


def np_savez(file_base, **kwargs):
    """
    reimplement the faulty (or so it appears) np.savez function
    """
    for k in kwargs:
        file_name = '%s.%s.gz' % (file_base, k)
        F = gzip.open(file_name, 'w')
        cPickle.dump(kwargs[k], F)
        F.close()

def np_load(file_base, what):
    """
    reimplement the crappy np.savez function
    """
    l.debug("loading %s %s" % (file_base, what))
    file_name = '%s.%s.gz' % (file_base, what)
    F = gzip.open(file_name, 'r')
    data = cPickle.load(F)
    F.close()
    return data


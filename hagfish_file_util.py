import os
import bz2
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

def np_exists(file_base, what):
    fn = '%s.%s.gz' % (file_base, what)
    return os.path.exists(fn)

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



def fastareader(f):
    if type(f) == type('hi'):
        #it is probably a filename
        if f[-4:] == '.bz2':
            F = bz2.BZ2File(f)
        else:
            F = open(f, 'r')
    else:
        #probably a file handle, or stdin
        F = f

    name, seq = "", []
    while True:
        l = F.readline()
        if not l: break
        
        l = l.strip()
        if not l: continue

        if l[0] == '>':
            if name and seq:
                yield name, "".join(seq)
            seq = []
            name = l[1:].split()[0]
        else:
            seq.append("".join(l.split()).lower())

    if name and seq:
        yield name, "".join(seq)

    F.close()

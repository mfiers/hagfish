import gzip
import cPickle

def np_savez(file_base, **kwargs):
    """
    reimplement the faulty (or so it appears) np.savez function
    """
    for k in kwargs:
        file_name = '%s.%s.gz' % (file_base, k)
        with gzip.open(file_name, 'w') as F:
            cPickle.dump(kwargs[k], F)

def np_load(file_base, what):
    """
    reimplement the crappy np.savez function
    """
    file_name = '%s.%s.gz' % (file_base, what)
    with gzip.open(file_name, 'r') as F:
        data = cPickle.load(F)
    return data


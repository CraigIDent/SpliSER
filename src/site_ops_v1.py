#site_ops.py
import pysam
from operator import add, truediv, mul, sub
from Gene_Site_Iter_Graph_v1 import Site



def trueDivCatchZero(array1, array2):
    """
    Takes an Array (array1) and divides it by another array (array2), elementwise.
    Where an element of array2 is equal to zero, the resulting value will be Zero (rather than undefined).
    """

    array3 = [0.0]*int(len(array1))
    for i, a2 in enumerate(array2):
        if a2 > 0.0:
            array3[i] = truediv(array1[i],a2)
    return array3


def calculateSSE(site):

    alpha = list(site.getAlphaCounts())
    beta1 = site.getBeta1Counts()
    beta2Simple = site.getBeta2SimpleCounts()
    betas = [x + y for x, y in zip(beta1, beta2Simple)]

    #if isbeta2Cryptic:
    #    beta2w = site.getBeta2WeightedCounts()
    #    betas = [x + y for x, y in zip(betas, beta2w)]

    denominator = [x + y for x, y in zip(alpha, betas)]

    site.setSSEs(trueDivCatchZero(list(alpha), list(denominator)))


def makeSingleSpliceSite(iChrom, iPos, numSamples, iStrand, isStranded):
    return(Site(
                chromosome = iChrom,
                pos = iPos,
                samples = numSamples,
                strand = iStrand,
                source = '',
                isStranded=isStranded
                )
            )









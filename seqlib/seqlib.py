#!/usr/bin/env python

"""
seqlib library for class assignment
"""

import copy
import numpy as np
import pandas as pd


class Seqlib:
    def __init__(self, ninds, nsites):  # no minfreq and maxfreq?
        self.ninds = ninds
        self.nsites = nsites
        # self.maxfreq = maxfreq    #this did not work...
        # self.minfreq = minfreq
        # self.arr = None
        # self.PD = None
        self.seqs = self._simulate()
        # originally had ninds and nsites inside () and no preceding _
        self.maf = self._get_maf()
        # did not even have anything like this... well, minfreq maybe?

    def _mutate(self, base):
        diff = set("ACTG") - set(base)
        return np.random.choice(list(diff))
        # made this a private fxn

    def _simulate(self):
        oseq = np.random.choice(list("ACGT"), size=self.nsites)
        arr = np.array([oseq for i in range(self.ninds)])
        muts = np.random.binomial(1, 0.1, (self.ninds, self.nsites))
        for col in range(self.nsites):
            newbase = self._mutate(arr[0, col])
            # made the above line self._mutate() instead of just mutate()
            mask = muts[:, col].astype(bool)
            arr[:, col][mask] = newbase
        missing = np.random.binomial(1, 0.1, (self.ninds, self.nsites))
        arr[missing.astype(bool)] = "N"
        return arr
        # self.arr = arr #this did not work.
        # using return arr means init will print out arr? or pass > next fxn?
        # I also made this fxn private


#    def maf(self):
#        arr = self.simulate()  #gen array
#        freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0] #adds up number of times the first seq matches the rest in the array
#        maf = freqs.copy()
#        maf[maf > 0.5] = 1 - maf[maf > 0.5]               #get minor allele freq. if freq cal'd above is >0.5, do 1-freq
#        return freqs

# not even sure what I was doing above...

    def _get_maf(self):
        maf = np.zeros(self.nsites)
        for col in range(self.nsites):
            thiscol = self.seqs[:, col]
            nmask = thiscol != "N"
            no_N_len = np.sum(nmask)
            first_non_n_base = thiscol[nmask][0]
            freq = np.sum(thiscol[nmask] != first_non_n_base) / no_N_len
            # the above line calculates frequency:
            # number of bases that deviate from first base in column
            # divided by # of bases in col, excluding Ns
            if freq > 0.5:
                maf[col] = 1 - freq
            else:
                maf[col] = freq
            # this if else statement assigns the minor allele freq only to maf
        return maf

        # i did not have anything like the above function in my package. lol

    def _filter_missing(self, maxmissing):
        freqmissing = np.sum(self.seqs == "N", axis=0) / self.seqs.shape[0]
        return freqmissing > maxmissing
        # made fxn rpviate + replaced two instances of "arr" with "self.seqs"
       
    def _filter_maf(self, minmaf):
        # freqs = np.sum(arr != arr[0], axis=0) / arr.shape[0]
        # maf = freqs.copy()
        # maf[maf > 0.5] = 1 - maf[maf > 0.5]
        # return arr[:, maf > minfreq]
        return self.maf < minmaf

# public fxns...
    def filter(self, minmaf, maxmissing):
        # maf = self.filter_maf(self.filter_missing(maxfreq), minfreq)
        # return maf
        """
        Applies maf and missing filters to the array 

        Parameters
        ----------
        minmaf: float
            The minimum minor allele frequency. Filter columns below this.
        maxmissing: float
            The maximum prop. missing data. Filter columns with prop Ns > this.
        """
        filter1 = self._filter_maf(minmaf)
        filter2 = self._filter_missing(maxmissing)
        fullfilter = filter1 + filter2
        return self.seqs[:, np.invert(fullfilter)]
        # apply both QC filters and return.

    def filter_seqlib(self, minmaf, maxmissing):
        """
        Applies maf and missing filters to the array and returns a copy
        of the seqlib object where the .seqs array has been filtered

        Parameters
        ----------
        minmaf: float
            The minimum minor allele frequency. Filter columns below this.
        maxmissing: float
            The maximum prop. missing data. Filter columns with prop Ns > this.
        """
        ## apply filters to get new array size
        newseqs = self.filter(minmaf, maxmissing)

        ## make a new copy of the seqlib object
        newself = copy.deepcopy(self)
        # deepcopy = copy embedded objects in self?
        newself.__init__(newseqs.shape[0], newseqs.shape[1])

        ## store the array (overwrite it)
        newself.seqs = newseqs

        ## call the _get_maf to match new array
        newself._get_maf()
        return newself

    def calculate_statistics(self):
        """ 
        Returns a dataframe of statistics on the seqs array. The earlier
        example from the notebook had a bug where var and inv were switched.
        """
        if self.seqs.size:
            nd = np.var(self.seqs == self.seqs[0], axis=0).mean()
            mf = np.mean(
                np.sum(self.seqs != self.seqs[0], axis=0) / self.seqs.shape[0])
            inv = np.all(self.seqs == self.seqs[0], axis=0).sum()
            var = self.seqs.shape[1] - inv
            return pd.Series(
                {"mean nucleotide diversity": nd,
                 "mean minor allele frequency": mf,
                 "invariant sites": inv,
                 "variable sites": var,
                })
        else:
            print("seqs array is empty")

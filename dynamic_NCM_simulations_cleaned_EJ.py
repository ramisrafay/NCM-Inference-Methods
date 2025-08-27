#!/usr/bin/env python3
#
# Companion code for the manuscript "How to quantify immigration from community
# abundance data using the Neutral Community Model", by Ramis Rafay, Eric
# Jones, David Sivak, and Jane Fowler.
#
# Code written by Eric Jones (jones.eric93@gmail.com)
#
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
import scipy
import csv
import pickle
import time
import bisect

################
# HELPER CLASSES
################

class Island():
    """State variable for an island population, using the island/mainland
    formulation of the Neutral Community Model. """
    def __init__(self, relative_abundances, Nt, m, singlet_initialization=False):
        """Initialize the state of the island based on input parameters.
        Inputs: self ... Island class variable
            relative_abundances ... relative abundances of species on mainland.
                                    length-S np.array w/ sum=1; S is # species.
            Nt ... integer, number of sites on the island
            m ... float between 0 and 1, immigration rate
            singlet_initialization ... boolean, initialize island with only the
                                       most abundant species
        """
        # initialize island-specific parameters
        self.pis = np.copy(relative_abundances)
        self.Nt = Nt
        self.m = m

        # default initialization method:
        # evenly allocate individuals based on mainland relative abundances 
        island_SAD = []
        for species,prob in enumerate(relative_abundances):
            island_SAD.append(int(np.round(self.Nt*prob)))
        # leftovers or leftunders deducted or added to most abundant species
        while sum(island_SAD) < self.Nt:
            island_SAD[0] += 1
        while sum(island_SAD) > self.Nt:
            island_SAD[0] -= 1
        # make sure we include rare species, even if none were initialized
        island_SAD = [int(xi) for xi in island_SAD]
        while len(island_SAD) < len(relative_abundances):
            island_SAD.append(0)
        # matching main text notation, xi is the island species abundance dist.
        self.xis = np.copy(island_SAD)

        # alternatively, initialize the island with only one species
        # (good for watching relaxation dynamics)
        if singlet_initialization:
            self.xis = [0 for i in range(len(self.xis))]
            self.xis[0] = Nt

    def evolve(self, num_steps=1, verbose=False):
        """Evolve the state of the island, following the update rule of the
        Neutral Community Model.
        Inputs: self ... Island class variable; the state to update
                num_steps ... integer, number of timesteps to update
                verbose ... boolean, whether to print log messages
        """

        t0 = time.time()

        # survey is a "flattened" SAD: the species of every observed individual
        # is enumerated
        survey = []
        for species,quantity in enumerate(self.xis):
            survey += [species]*quantity

        # preload the "random individuals" to be removed and introduced
        deletions = np.random.choice(sum(self.xis), size=num_steps)
        additions = np.random.choice(sum(self.xis)-1, size=num_steps)
        immigrants = np.random.choice(len(self.pis), p=self.pis, size=num_steps)
        immigration_draws = np.random.random(size=num_steps)

        # NCM update timestep
        for ii in range(num_steps):
            if verbose:
                self.xis = np.bincount(survey)
                print(self.xis)
                print()
            removed = survey.pop(deletions[ii])
            if verbose: print('removing species {}'.format(removed))

            immigration_draw = immigration_draws[ii]
            if immigration_draw < self.m:
                insert_idx = bisect.bisect_left(survey, immigrants[ii])
                survey.insert(insert_idx, immigrants[ii])
                if verbose: print('immigration: adding species {}'.format(
                                  immigrants[ii]))
            else:
                survey.insert(additions[ii], survey[additions[ii]])
                if verbose: print('on island: adding species {}'.format(
                                  additions[ii]))

        self.xis = np.bincount(survey)
        self.xis = [int(xi) for xi in self.xis]
        # make sure rare species (with 0 abundance) are included
        while len(self.xis) < len(self.pis):
            self.xis.append(0)

        if verbose:
            print(self.xis)
            print('t={}s'.format(time.time() - t0))
            print()

    def evolve_allee(self, num_steps=1, verbose=False, allee_threshold=1):
        """Evolve the state of the island, based on the update rule of the
        Neutral Community model but include an Allee effect. Most code is the
        same as in evolve(). Note evolve_allee() takes much longer than
        evolve() to run, expecially with moderate allee_threshold.
        Inputs: self ... Island class variable; the state to update
                num_steps ... integer, number of timesteps to update
                verbose ... boolean, whether to print log messages
        """
        t0 = time.time()
        survey = []
        for species,quantity in enumerate(self.xis):
            survey += [species]*quantity

        deletions = np.random.choice(sum(self.xis), size=num_steps)
        additions = np.random.choice(sum(self.xis)-1, size=num_steps)
        immigrants = np.random.choice(len(self.pis), p=self.pis, size=num_steps)
        immigration_draws = np.random.random(size=num_steps)
        N = len(survey)

        for ii in range(num_steps):
            if verbose: print(self.xis)
            removed = survey.pop(deletions[ii])
            if verbose: print('removing species {}'.format(removed))

            immigration_draw = immigration_draws[ii]
            if immigration_draw < self.m: # replaced by immigration
                # each draw, calculate the drawn species'  abundance
                lb = bisect.bisect_left(survey, immigrants[ii])
                ub = bisect.bisect_right(survey, immigrants[ii])
                my_abundance = ub - lb
                # Allee effect step:
                if my_abundance >= allee_threshold: # colonization allowed
                    insert_idx = bisect.bisect_left(survey, immigrants[ii])
                    survey.insert(insert_idx, immigrants[ii])
                    if verbose: print('imm: 1: added species {} w abundance {}'
                                      .format( immigrants[ii], my_abundance))
                else: # colonization forbidden, draw a new "allowed" species
                    if verbose: print('imm: fail to add species {} abundance {}'
                                      .format(immigrants[ii], my_abundance))
                    # construct a list of species that are "allowed" to be
                    # drawn, i.e. have abundance >= allee_threshold
                    temp_xis = np.bincount(survey)
                    temp_survey = []
                    temp_pis = []
                    for species,quantity in enumerate(temp_xis):
                        if quantity >= allee_threshold:
                            temp_pis.append(self.pis[species])
                        else:
                            temp_pis.append(0)
                    temp_pis = [val/sum(temp_pis) for val in temp_pis]
                    # select and insert one of these "allowed" species
                    addition = np.random.choice(len(temp_pis), p=temp_pis)
                    insert_idx = bisect.bisect_left(survey, addition)
                    survey.insert(insert_idx, addition)
                    if verbose: print('imm: 2: added species {} w pi {}'
                                      .format(addition, temp_pis[addition]))
            else: # replaced by species on island
                # calculate drawn species' abundance
                lb = bisect.bisect_left(survey, survey[additions[ii]])
                ub = bisect.bisect_right(survey, survey[additions[ii]])
                my_abundance = ub - lb
                # Allee effect step:
                if my_abundance >= allee_threshold:
                    survey.insert(additions[ii], survey[additions[ii]])
                    if verbose: print('island: 1: added species {} w abundance {}'
                                      .format(survey[immigrants[ii]], my_abundance))
                else:
                    if verbose: print('island: fail add species {} abundance {}'
                                      .format(survey[additions[ii]], my_abundance))
                    temp_xis = np.bincount(survey)
                    temp_survey = []
                    for species,quantity in enumerate(temp_xis):
                        if quantity >= allee_threshold:
                            temp_survey += [species]*quantity
                    addition = np.random.choice(temp_survey)
                    addition_idx = bisect.bisect_left(survey, addition)
                    survey.insert(addition_idx, addition)
                    if verbose: print('island: 2: added species {} w abundance {}'
                                      .format(addition, temp_xis[addition]))

        self.xis = np.bincount(survey)
        self.xis = [int(xi) for xi in self.xis]
        while len(self.xis) < len(self.pis):
            self.xis.append(0)

        if verbose:
            print(self.xis)
            print('t={}s'.format(time.time() - t0))
            print()

class MetaIsland():
    """State variable for an metapopulation of many islands. Variant of the
    island/mainland Neutral Community Model where when one island updates, all
    other islands act as the mainland. Code otherwise similar to the Island
    class.  """
    def __init__(self, relative_abundances, Nt, m, num_islands):
        """Initialize the state of the metapopulation based on input parameters.
        Inputs: self ... Island class variable
            relative_abundances ... relative abundances of species on mainland.
                                    length-S np.array w/ sum=1; S is # species.
            Nt ... integer, number of sites on each island
            m ... float between 0 and 1, immigration rate
            num_islands ... number of islands that compose the metapopulation
        """

        self.relative_abundances = relative_abundances
        self.Nt = Nt
        self.m = m
        self.num_islands = num_islands
        self.islands = {j: [] for j in range(self.num_islands)}

        island_SAD = []
        for species,prob in enumerate(self.relative_abundances):
            island_SAD.append(int(np.round(self.Nt*prob)))
        while sum(island_SAD) < self.Nt:
            island_SAD[0] += 1
        while sum(island_SAD) > self.Nt:
            island_SAD[0] -= 1
        island_SAD = [int(xi) for xi in island_SAD]
        while len(island_SAD) < len(relative_abundances):
            island_SAD.append(0)

        # islands is a dictionary of SADs
        for island in range(num_islands):
            self.islands[island] = [int(xi) for xi in island_SAD]

    def evolve(self, num_steps=1, verbose=False):
        """Evolve the state of the metapopulation, following the modified NCM
        rules for metapopulations. 
        Inputs: self ... Island class variable; the state to update
                num_steps ... integer, number of timesteps to update
                verbose ... boolean, whether to print log messages
        """
        t0 = time.time()
        survey = [] # survey is a list of all individuals across all islands
        for island in range(self.num_islands):
            for species,quantity in enumerate(self.islands[island]):
                survey += [species]*quantity

        # proload random numbers needed for update steps
        # which individuals will be removed each timestep:
        deletions = np.random.choice(len(survey), size=num_steps)
        # the islands on which these removed individuals are located:
        focal_islands = [int(val/self.Nt) for val in deletions]
        # if replaced by individuals from focal island:
        additions = np.random.choice(self.Nt-1, size=num_steps)
        # if replaced by individuals from random other island:
        drawn_islands = np.random.choice(range(1, self.num_islands),size=num_steps)
        # the immigrant from this drawn_island:
        immigrants = np.random.choice(self.Nt, size=num_steps)
        immigration_draws = np.random.random(size=num_steps)

        # update step: all modifications happen in 'survey' variable
        for ii in range(num_steps):
            removed = survey.pop(deletions[ii])
            focal_island = focal_islands[ii]
            if verbose: print('removing species {} from island {}'.format(
                              removed, focal_island))

            immigration_draw = immigration_draws[ii]
            if immigration_draw < self.m:
                drawn_island = (focal_island + drawn_islands[ii]) % self.num_islands
                addition = survey[self.Nt*drawn_island + additions[ii]]
                insert_idx = (bisect.bisect_left(
                    survey[self.Nt*focal_island:self.Nt*(focal_island+1)-1],
                    addition))
                survey.insert(insert_idx + self.Nt*focal_island, addition)
                if verbose: print(
                        'on island {} (immigration): adding species {} from island {}'.format(
                            focal_island, addition, drawn_island))
            else:
                survey.insert(self.Nt*focal_island + additions[ii],
                              survey[self.Nt*focal_island + additions[ii]])
                if verbose: print('on island {} (self): adding species {}'.format(
                                  focal_island,
                                  survey[self.Nt*focal_island + additions[ii]]))

        # calculate SADs from each island
        for island in range(self.num_islands):
            temp_isl = np.bincount(survey[self.Nt*island:self.Nt*(island+1)])
            temp_isl = [int(xi) for xi in temp_isl]
            while len(temp_isl) < len(self.relative_abundances):
                temp_isl.append(0)
            self.islands[island] = temp_isl

        if verbose:
            print(self.islands)
            print('t={}s'.format(time.time() - t0))
            print()

##################
# HELPER FUNCTIONS 
##################

def import_species_abundance_distribution():
    """Import the Barro Colorado Island abundance table, and use it as mainland
    relative abundances.
    Outputs: relative_abundances ... relative abundances of species on mainland.
                                     length-S np.array w/ sum=1; S is # species.
        all_reads: J x S np.array; total abundance table for BCI;
                   each row is a sample, each column is a species
    """
    all_reads = []
    with open('./condit_2002_abundance_table.csv') as csvfile:
        spamreader = csv.reader(csvfile)
        for i,row in enumerate(spamreader):
            if i == 0: continue
            row_reads = []
            for j,elem in enumerate(row):
                if j == 0: continue
                row_reads.append(float(elem))
            all_reads.append(row_reads)

    all_reads = np.array(all_reads).T
    species_counts = np.sum(all_reads, axis=0)
    total_reads = np.sum(species_counts)
    relative_abundances = species_counts/total_reads
    relative_abundances = [val for val in relative_abundances if val > 0]
    relative_abundances = np.array(sorted(relative_abundances)[::-1])

    return relative_abundances, all_reads

def get_SAD(relative_abundances, Nt, m, J, t_end, R=0, keyword='',
            allee_threshold=None):
    """Load and/or generate a Species Abundance Distribution for an island
    obeying the NCM.
    Inputs: 
        relative_abundances ... relative abundances of species on mainland.
                                length-S np.array w/ sum=1; S is # species.
        Nt ... integer, # sites on island
        m ... float between 0 and 1, immigration rate
        J ... integer, # simulations to perform (i.e., # events in histogram)
        t_end ... integer, number of timesteps to simulate
        R ... integer, read depth
        keyword ... float, convey additional attributes
        allee_threshold ... int, threshold for Allee effect population size
    Outputs:
        SAD ... J x S array, rows are different island simulations,
                columns are different species
    Usage: SAD[2][5] = abundance of 6th-most-abundant species in simulation 2
    """

    filename = 'vars/SAD_{}_Nt_{}_t_end_{}_m_{}_J_{}_R_{}.pi'.format(
                keyword, Nt, t_end, m, J, R)
    if keyword == 'allee':
        filename = 'vars/SAD_{}_Nt_{}_t_end_{}_m_{}_J_{}_R_{}_allee_{}.pi'.format(
                    keyword, Nt, t_end, m, J, R, allee_threshold)
    try:
        with open(filename, 'rb') as f:
            print('loading {} || Nt={}, m={}, J={}, t_end={}, R={}'.format(
                keyword, Nt, m, J, t_end, R))
            SAD = pickle.load(f)
    except FileNotFoundError:
        print('simulating {} || Nt={}, m={}, J={}, t_end={}, R={}'.format(
            keyword, Nt, m, J, t_end, R))
        SAD = []
        for i in range(J):
            if i % 20 == 0: print('... {}'.format(i))
            if keyword == 'singlet_no_sampling':
                isl = Island(relative_abundances, Nt, m,
                             singlet_initialization=True)
            else:
                isl = Island(relative_abundances, Nt, m,
                             singlet_initialization=False)

            if keyword == 'allee':
                isl.evolve_allee(num_steps=t_end,
                                 allee_threshold=allee_threshold,
                                 verbose=False)
            else:
                isl.evolve(num_steps=t_end)
            SAD.append(isl.xis)
        with open(filename, 'wb') as f:
            pickle.dump(SAD, f)
    return SAD

def get_meta_SAD(relative_abundances, Nt, m, J, num_islands, t_end, R=0):
    """Load and/or generate Species Abundance Distributions for a
    metacommunity (i.e., num_islands islands each obeying the NCM).
    Inputs: 
        relative_abundances ... relative abundances of species on mainland.
                                length-S np.array w/ sum=1; S is # species.
        Nt ... integer, # sites on island
        m ... float between 0 and 1, immigration rate
        J ... integer, # simulations to perform (i.e., # events in histogram)
        num_islands ... integer, number of islands in the metacommunity
        t_end ... integer, number of timesteps to simulate
        R ... integer, read depth
        keyword ... float, convey additional attributes
        allee_threshold ... int, threshold for Allee effect population size
    Outputs:
        SAD ... dict of SADs, keys are 0 ... num_islands
    Usage: SAD[13][2][5] = abundance on island 13 of 6th-most-abundant species
                           in simulation 2
    """
    keyword = 'meta'
    filename = 'vars/SAD_{}_Nt_{}_t_end_{}_m_{}_J_{}_isl_{}_R_{}.pi'.format(
                keyword, Nt, t_end, m, J, num_islands, R)
    try:
        with open(filename, 'rb') as f:
            print('loading {} || Nt={}, m={}, J={}, num isl={}, t_end={}, R={}'.format(
                keyword, Nt, m, J, num_islands, t_end, R))
            SAD = pickle.load(f)
    except FileNotFoundError:
        print('simulating {} || Nt={}, m={}, J={}, num isl={}, t_end={}, R={}'.format(
            keyword, Nt, m, J, num_islands, t_end, R))
        SAD = {island: [] for island in range(num_islands)}
        for i in range(J):
            if i % 100 == 0: print('... {}'.format(i))
            isl = MetaIsland(relative_abundances, Nt, m, num_islands)
            isl.evolve(num_steps=t_end)

            for island in range(num_islands):
                SAD[island].append(isl.islands[island])
        with open(filename, 'wb') as f:
            pickle.dump(SAD, f)
    return SAD

################################
# GENERATE SUPPLEMENTARY FIGURES
################################

def make_figure_S12(verbose=False):
    """Generates Figure S12 of the manuscript by running NCM simulations with
    an Allee effect, then inferring Ntm using the three inference methods as a
    function of read depth."""

    relative_abundances, all_reads = import_species_abundance_distribution()
    Nt = 25000
    m = 0.02 # hence Nt*m = 500
    J = 500 # num samples
    J_subsample = 50 # the "sample size S" of a single sampling effort
                     # hence we bootstrap over 10 collection efforts.
    t_end = 2e6; t_end = int(t_end)
    allee_list = [0, 1, 5, 10] # list of Allee effects to simulate
    R_list = np.logspace(np.log10(100), np.log10(10000), 20)
    R_list = [int(val) for val in R_list] # list of read depths for sampling

    colors = plt.cm.jet(np.linspace(0, 1, len(allee_list)))
    filename = 'vars/bs_infer_Rs_{}_Rm_{}_Nt_{}_J_{}_Jsub_{}_t_end_{}.pi'.format(
                len(R_list), min(R_list), Nt, J, J_subsample, t_end)

    # load Allee effect simulations, if possible
    try:
        with open(filename, 'rb') as f:
            if verbose: print('loading {}'.format(filename))
            allee_list, R_list, Ntms_inferred = pickle.load(f)
            # Usage:
            # Ntms_occ, Ntms_vars, Ntms_dmll = NTMs_inferred
            # np.shape(Ntms_occ) = (len(allee_thresholds), len(R_lists))
    except FileNotFoundError:
        # initialize inferred Ntm values by method:
        Ntms_occ = []
        Ntms_var = []
        Ntms_dmll = []
        for jj,allee_threshold in enumerate(allee_list[3:]):
            Ntm_occ_allee = []
            Ntm_var_allee = []
            Ntm_dmll_allee = []
            master_SAD = get_SAD(relative_abundances, Nt, m, J, t_end,
                          keyword='allee', allee_threshold=allee_threshold)
            # bootstrapping step: iterating over J_subsample-sized loops
            for SAD_idx in range(int(len(master_SAD)/J_subsample)):
                # SAD for a single subsample:
                SAD = master_SAD[SAD_idx*J_subsample:(SAD_idx+1)*J_subsample]
                SAD_main = np.copy(SAD)
                Ntm_occ = []
                Ntm_var = []
                Ntm_dmll = []
                for R in R_list: # iterating over read depths
                    if verbose: print('allee threshold = {}, SAD_idx = {}, R = {}'.format(
                          allee_threshold, SAD_idx, R))
                    new_SAD = []
                    # sampling from the full island of size Nt. R << Nt.
                    for sample in SAD_main:
                        probs = [val/sum(sample) for val in sample]
                        #draw = np.random.multinomial(sum(sample), probs)
                        draw = np.random.multinomial(R, probs)
                        new_SAD.append(draw)
                    SAD = new_SAD

                    xis = [[val/sum(SAD[0]) for val in sample] for sample in SAD]
                    yis = [[val for val in sample] for sample in SAD]
                    pis = relative_abundances 

                    def occupancy_method(Ntm):
                        """Returns the error for the occupancy-based inference
                        method. Relies on having xi and pi variables in scope.
                        Inputs: Ntm ... float, input from optimization function
                        Outputs: utility ... float, quantity to minimize
                        """
                        d = 1/R # limit of detection
                        # calculating empirical and theoretical occupancies  
                        emp_occs = np.mean([[1 if val else 0 for val in sample] for sample
                                            in xis], axis=0)
                        theor_occs = [scipy.special.betaincc(Ntm*pi, Ntm*(1-pi), d)
                                      for pi in pis]
                        theor_occs = np.array([val if not np.isnan(val) else 0 for val in
                                      theor_occs])
                        res = theor_occs - emp_occs
                        utility = np.sqrt(sum([xx*xx for xx in res]))
                        return utility 

                    # infer Ntm with the occupancy method
                    res = scipy.optimize.minimize_scalar(occupancy_method, bounds=[1, 1e6])
                    Ntm_occ.append(res.x)

                    def var_method(Ntm):
                        """Returns the error for the variance-based inference
                        method. Relies on having R, xi, and pi variables in scope.
                        Inputs: Ntm ... float, input from optimization function
                        Outputs: utility ... float, quantity to minimize
                        """
                        a = np.var(xis, axis=0)
                        b = np.array([pi*(1-pi)/(Ntm+1) + pi*(1-pi)/R for pi in pis])
                        res =  a - b
                        utility = np.sqrt(sum([xx*xx for xx in res]))
                        return utility 

                    # infer Ntm with the variance method
                    res = scipy.optimize.minimize_scalar(var_method)
                    Ntm_var.append(res.x)

                    def DMLL_method(Ntm):
                        """Returns the error for the DM-LL inference
                        method. Relies on having R, yi, and pi variables in scope.
                        Inputs: Ntm ... float, input from optimization function
                        Outputs: utility ... float, quantity to minimize
                        """
                        constant = (np.log(R) + scipy.special.loggamma(Ntm)
                                    + scipy.special.loggamma(R) 
                                    - scipy.special.loggamma(Ntm + R))
                        likelihood = 0
                        for sample in SAD:
                            summand = constant
                            for yi,pi in zip(sample, pis):
                                if yi == 0: continue
                                summand -= (np.log(yi) + scipy.special.loggamma(Ntm*pi)
                                            + scipy.special.loggamma(yi)
                                            - scipy.special.loggamma(Ntm*pi + yi))
                            likelihood += summand
                        return -likelihood

                    # infer Ntm with DM-LL
                    res = scipy.optimize.minimize_scalar(DMLL_method, bounds=[1, 1e7])
                    Ntm_dmll.append(res.x)
                    print('occupancy method: Ntm = {}'.format(Ntm_occ[-1]))
                    print('variance method: Ntm = {}'.format(Ntm_var[-1]))
                    print('DMLL method: Ntm = {}'.format(Ntm_dmll[-1]))
                    print()

                Ntm_var_allee.append(Ntm_var)
                Ntm_occ_allee.append(Ntm_occ)
                Ntm_dmll_allee.append(Ntm_dmll)

            Ntms_var.append(Ntm_var_allee)
            Ntms_occ.append(Ntm_occ_allee)
            Ntms_dmll.append(Ntm_dmll_allee)

        Ntms_inferred = [Ntms_occ, Ntms_var, Ntms_dmll]
        with open(filename, 'wb') as f:
            pickle.dump([allee_list, R_list, Ntms_inferred], f)

    # plot Fig S12:
    fig, ax = plt.subplots(ncols=3, figsize=(9.6, 4.8), sharey=True)
    colors = plt.cm.jet(np.linspace(0, 1, len(allee_list)))
    inference_types = ['occupancy', 'variance', 'DM-LL']
    for jj,Ntms in enumerate(Ntms_inferred):
        for ii,reps in enumerate(Ntms):
            means = []; lbs = []; ubs = []
            reps = np.array(reps).T
            R_vals = []
            for R_val,sample in zip(R_list, reps):
                # QC for inferred Ntm values
                sample = [val for val in sample if val < 10000]
                if len(sample) == 0:
                    continue

                R_vals.append(R_val)
                mean = np.mean(sample)
                val = scipy.stats.bootstrap((sample,), np.mean)
                lb = val.confidence_interval[0]
                ub = val.confidence_interval[1]
                means.append(mean)
                if (mean > 1000) or (mean < 1.01):
                    lbs.append(mean)
                    ubs.append(mean)
                else:
                    lbs.append(lb)
                    ubs.append(ub)

            ax[jj].errorbar(R_vals, means, fmt='o',
                            yerr=[[mean - lb for mean,lb in zip(means, lbs)],
                                  [ub - mean for mean,ub in zip(means, ubs)]],
                            color=colors[ii], ms=7, label=allee_list[ii])
            ax[jj].set_title(inference_types[jj])
        ax[jj].axhline(y=500, color='k')
        ax[jj].axis([100, 10000, 0, 1000])
        ax[jj].set_xscale('log')
        ax[jj].set_xlabel(r'read depth $R$')
    ax[0].set_ylabel(r'inferred $N_t m$')
    ax[0].legend(title='Allee threshold', fontsize=10, title_fontsize=10)
    plt.savefig('figs/Fig_S12_Nt_{}_Nmax_{}.pdf'.format(
        Nt, max(R_list)), bbox_inches='tight')

def make_figure_S1(verbose=False):
    """Generate manuscript Figure S1."""
    relative_abundances, all_reads = import_species_abundance_distribution()
    Nt = 500; num_islands = 50;
    J = 5000 # J is number of simulations run
    species_rank = 1 # which species' abundance distribution should be plotted?

    
    fig, [ax1, ax2, ax3] = plt.subplots(ncols=3, figsize=(12.8, 4.8))

    # ax1: Fig S1a: Plotting the ``mainland'' rank abundance distribution
    if verbose: print('Generating Fig S1a')
    ax1.plot(range(1, len(relative_abundances) + 1), relative_abundances, 'k.')
    ax1.set_xlabel(r'rank $i$')
    ax1.set_ylabel(r'relative abundance $p_i$')
    ax1.set_yscale('log')
    ax1.set_xticks([1, 100, 200])
    ax1.text(-80, 0.1, r'{\bf \textsf{a})}', fontsize=20)

    # ax2: Fig S1b: Dynamic NCM simulations match Beta distribution theory
    if verbose: print('Generating Fig S1b')
    m = 3*0.1
    bins = np.linspace(0, 0.2, 100)
    t_end = int(5e5)
    SAD = get_SAD(relative_abundances, Nt, m, J, t_end,
                  keyword='singlet_no_sampling')
    yis = [[val/sum(SAD[0]) for val in sample] for sample in SAD]
    rel_abundances = [val[species_rank] for val in yis]
    pi = np.mean(rel_abundances)
    # plot histograms from dynamic simulations
    ax2.hist(rel_abundances, bins=bins, histtype='step',
            label='dynamic NCM simulations', lw=3, density=True,
            color='grey')
    # plot theoretical beta distribution
    xs = np.linspace(0, 0.2, 201)
    alpha = Nt*m*pi; beta = Nt*m*(1 - pi)
    ys = []
    for x in xs:
        try:
            log_ans = (-scipy.special.betaln(alpha, beta) +
                       (alpha-1)*np.log(x) + (beta-1)*np.log(1-x))
            ys.append(np.exp(log_ans))
        except RuntimeWarning:
            ys.append(np.nan)
    ax2.plot(xs, ys, 'k', lw=3,
            label=r'Beta($N_{\scriptsize\textsf{T}} m p_i,\ N_{\scriptsize\textsf{T}} m (1 - p_i)$)',
            zorder=10)

    # ax3: Fig S1c: Metapopulation NCM simulations match Beta distribution
    if verbose: print('Generating Fig S1c')
    t_ends = [0, 1e3, 5e3, 1e4, 5e4, 1e5] # to watch SAD relaxation
    num_t_ends = len(t_ends)
    m = 3*0.1
    bins = np.linspace(0, 0.2, 100)
    colors = plt.cm.jet(np.linspace(0, 1, num_t_ends))

    for ii,t_end in enumerate(t_ends):
        t_end = int(t_end)
        SAD = get_meta_SAD(relative_abundances, Nt, m, J, num_islands, t_end)
        # plot histograms from dynamic simulations
        for jj,island in enumerate(range(num_islands)):
            if verbose and jj % 10 == 0:
                print(f'surveying island {jj} for t_end {t_end}')
            yis = [[val/sum(SAD[island][0]) for val in sample]
                   for sample in SAD[island]]
            rel_abundances = [val[species_rank] for val in yis]
            if jj == 0:
                ax3.hist(rel_abundances, bins=bins, histtype='step',
                         label='{:,}'.format(t_end), lw=0.8, density=True,
                        color=colors[ii])
            else:
                ax3.hist(rel_abundances, bins=bins, histtype='step',
                        lw=0.8, density=True, color=colors[ii])
        # plot theoretical beta distribution
        if t_end == max(t_ends):
            pi = relative_abundances[species_rank]
            xs = np.linspace(0, 0.2, 201)
            alpha = Nt*m*pi; beta = Nt*m*(1 - pi)
            ys = []
            for x in xs:
                try:
                    log_ans = (-scipy.special.betaln(alpha, beta) +
                               (alpha-1)*np.log(x) + (beta-1)*np.log(1-x))
                    ys.append(np.exp(log_ans))
                except RuntimeWarning:
                    ys.append(np.nan)
            ax3.plot(xs, ys, 'k', lw=3, zorder=20, label='Beta dist.')

    ax2.axis([0, 0.15, 0, 40])
    ax3.axis([0, 0.15, 0, 40])
    ax2.set_xlabel('relative abundance')
    ax3.set_xlabel('relative abundance')
    ax2.set_ylabel('probability')
    ax3.set_ylabel('probability')
    ax2.set_xticks([0, 0.05, 0.1, 0.15])
    ax3.set_xticks([0, 0.05, 0.1, 0.15])
    ax2.text(-0.029, 39.3, r'{\bf \textsf{b})}', fontsize=20)
    ax3.text(-0.029, 39.3, r'{\bf \textsf{c})}', fontsize=20)
    ax2.legend(fontsize=11)
    ax3.legend(title='elapsed timesteps', title_fontsize=11, fontsize=11)

    plt.tight_layout(pad=0.5)
    plt.savefig('figs/Fig_S1_Nt_{}_m_{}_J_{}_isl_{}_rank_{}.pdf'.format(
                Nt, m, J, num_islands, species_rank), bbox_inches='tight')


################################
# USER INPUT FOR NCM SIMULATIONS
################################

if __name__ == "__main__":
    # MINIMAL WORKING EXAMPLE FOR SIMULATING NCM
    relative_abundances = [0.3, 0.3, 0.1, 0.1, 0.1, 0.1]
    Nt = 500
    m = 0.5
    isl = Island(relative_abundances, Nt, m)
    isl.evolve(num_steps=5, verbose=True)

    # GENERATE FIGURES S1 AND S12
    # make_figure_S1(verbose=True) 
    # make_figure_S12(verbose=True)
    # NOTE: these functions take 5ish hours to run; for speedup lower J

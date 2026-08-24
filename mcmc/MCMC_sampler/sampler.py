# This code has been taken from CosmoHammer (proper citation to be made in the paper when published)
# with lots of Modification as required

import emcee
from multiprocessing import Pool
import numpy as np
import time
import logging

import file_configure as c
from utils import getLogger


class MCMCsampler():

    def __init__(self, params, likelihoodComputation, filePrefix, chain_storage_path, fileroot,
                 walkersRatio, sampleIterations, n_cores=1, backend = None, logLevel=logging.INFO):
        """
        MCMC sampler implementation.

        Parameters
        ----------
        params                 : array — shape (N, 4), columns: [value, min, max, width]
        likelihoodComputation  : callable — log-likelihood function
        filePrefix             : str — prefix for log file
        chain_storage_path     : str — directory to save chain files
        fileroot               : str — root name for chain files
        walkersRatio           : int — nwalkers = paramCount * walkersRatio
        sampleIterations       : int — number of MCMC steps
        n_cores                : int — number of CPU cores for parallelism (default 1 = serial)
        """
        self.params                      = params
        self.paramValues                 = self.params[:, 0]
        self.paramWidths                 = self.params[:, 3]
        self.likelihoodComputationChain  = likelihoodComputation
        self.walkersRatio                = walkersRatio
        self.filePrefix                  = filePrefix
        self.chain_storage_path          = chain_storage_path
        self.file_root                   = fileroot
        self.n_cores                     = n_cores
        self.paramCount                  = len(self.paramValues)
        self.nwalkers                    = self.paramCount * walkersRatio
        self.sampleIterations            = sampleIterations
        self.backend                     = backend

        if not hasattr(self.likelihoodComputationChain, "params"):
            self.likelihoodComputationChain.params = params

        # Set up logging
        self._configureLogging(filePrefix + c.LOG_FILE_SUFFIX, logLevel)

    # ------------------------------------------------------------------

    def InitPosGenerator(self):
        return [
            self.paramValues + np.random.normal(size=self.paramCount) * self.paramWidths
            for i in range(self.nwalkers)
        ]

    # ------------------------------------------------------------------

    def startSampling(self):
        for k in range(self.nwalkers):
            f = open(self.chain_storage_path + '/' + self.file_root + "_" + str(k + 1) + ".txt", "w+")
            f.write('# Weight, Lnprob, different params\n')
            f.close()

        pos = self.InitPosGenerator()

        self.log("start sampling — no burn-in. Remove burn-in steps when producing corner plots.")
        self.log(f"Using {self.n_cores} core(s).")
        start = time.time()

        if self.n_cores > 1:
            with Pool(processes=self.n_cores) as pool:
                sampler = self.createEmceeSampler(pool=pool)
                self.sample(sampler, pos)
        else:
            sampler = self.createEmceeSampler(pool=None)
            self.sample(sampler, pos)

        end = time.time()
        self.log("sampling done! Took: " + str(round(end - start, 4)) + "s")

    # ------------------------------------------------------------------

    def createEmceeSampler(self, pool=None):
        """
        Create the emcee sampler.

        pool=None  → serial execution
        pool=Pool  → parallel execution via multiprocessing
        """
        if self.backend is not None:

            return emcee.EnsembleSampler(
            self.nwalkers,
            self.paramCount,
            self.likelihoodComputationChain,
            pool=pool, 
            backend=self.backend)

        else:
            return emcee.EnsembleSampler(
            self.nwalkers,
            self.paramCount,
            self.likelihoodComputationChain,
            pool=pool)

    # ------------------------------------------------------------------
    
    '''
    def sample(self, sampler, InitPos):
        """
        Run the MCMC sampler and save chain to disk every save_steps iterations.
        """
        print('MCMC sampling with nwalkers=%s, Iterations=%s, Number of parameters=%s\n' % (
            self.nwalkers, self.sampleIterations, self.paramCount))

        nblobs         = 2
        samplenumber   = self.sampleIterations
        chain_arr      = np.empty([self.nwalkers, samplenumber, self.paramCount])
        lnprob_arr     = -np.inf * np.ones([self.nwalkers, self.sampleIterations])
        blobs_arr      = np.empty([self.nwalkers, self.sampleIterations, nblobs])

        step_index       = 0
        save_steps_start = 0
        counter          = 1
        save_steps       = 10

        for pos, prob, stat, blobs in sampler.sample(InitPos, iterations = self.sampleIterations):
            chain_arr[:, step_index, :] = pos
            lnprob_arr[:, step_index]   = prob
            blobs_arr[:, step_index]    = np.asarray(blobs)
            step_index += 1

            if np.remainder(step_index, save_steps) == 0:
                for k in range(self.nwalkers):
                    print('step, save_step', step_index, save_steps)
                    f = open(self.chain_storage_path + '/' + self.file_root + "_" + str(k + 1) + ".txt", "a")
                    for i in range(save_steps_start, save_steps_start + save_steps):
                        s  = "{0:6d}".format(1)   # dummy weight
                        s += " " + "{:.6e}".format(-lnprob_arr[k, i])
                        for kk in range(self.paramCount):
                            s += " " + "{:.6e}".format(chain_arr[k, i, kk])
                        for kk in range(len(blobs_arr[k, i, :])):
                            s += " " + "{:.6e}".format(blobs_arr[k, i, kk])
                        s += "\n"
                        f.write(s)
                    f.close()
                save_steps_start += save_steps

            if counter % 1 == 0:
                self.log("Iteration finished with total sample Number " + str(counter * self.nwalkers) + '\n')

            counter += 1
    '''
    
    def sample(self, sampler, InitPos):
        completed = self.backend.iteration

        if completed == 0:
            self.backend.reset(self.nwalkers, self.paramCount)
            start_state = InitPos
        else:
            self.log(f"Completed step before this run {completed}")
            print(f"Completed step before this run {completed}")
            start_state = self.backend.get_last_sample()      # Resume from backend

        
        remaining = self.sampleIterations - completed

        for state in sampler.sample(start_state,
                            iterations=remaining):

            if sampler.iteration % 10 == 0:
                self.log(
                         f"Completed {sampler.iteration}/{self.sampleIterations} iterations")

    # ------------------------------------------------------------------

    def log(self, message, level=logging.INFO):
        getLogger().log(level, message)

    # ------------------------------------------------------------------

    def _configureLogging(self, filename, logLevel):
        logger = getLogger()
        logger.setLevel(logLevel)
        fh = logging.FileHandler(filename, "w")
        fh.setLevel(logLevel)
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(asctime)s %(levelname)s:%(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
        logger.addHandler(fh)
        logger.addHandler(ch)

import os, sys, shutil
import numpy as np
import pandas as pd
from multiprocessing import  Pool
import subprocess
from LatticeData import *

'''
Notes and plans for changes

For even L=10 lattices, catalog_sites is extremely long. 
    Is there a way to submit a batch job to do this in parallel and then save the results? 
    Try parallelizing: https://dmnfarrell.github.io/python/parallelize-python-df
    What about finding a way to output this data in the simulation itself?
'''

class TestingSuite(LatticeData):
    
    def setup(self,run_dirpath, use_full_filepath = False, test_log = "test.log"):
        self.set_run_dir(run_dirpath)
        self.use_full_filepath = use_full_filepath
        self.set_test_logfile(test_log = test_log)
        self.site_catalog = None
        self.test_params = None
    
    def set_run_dir(self,run_dirpath):
        #select location of data
        if self.use_full_filepath:
            self.test_run_dir = run_dirpath+'/'
        else:
            self.test_run_dir = os.getcwd()+"/"+run_dirpath
        
        
    def set_test_logfile(self,test_log = "test.log"):
        self.test_log = test_log
        
    def find_site(self, i = 0, j = 0, step = 0):
        sites = self.site_catalog.copy()
        mask = (sites["step"] == step) & (sites["i"] == i) & (sites["j"] == j)
        return sites[mask]
                    
    def create_catalog(self, logfile_subset):
        
    def catalog_sites(self, **kwargs):
        missing_params = []
        for p in self.parameters:
            if p not in kwargs.keys():
                missing_params.append(p)
        if len(missing_params) > 0:
            self._error_message("Missing parameters in input: ", missing_params)
            return None
        self.test_params = kwargs
        len_file, testing_logfile = self._get_test_logfile()
        catalog = pd.DataFrame()#dict()
        start = self._search(testing_logfile,"Initializing lattice")
        MCstep = None
        i = None
        j = None
        for n,line in enumerate(testing_logfile[start:]):
            temp = line.strip('\n')
            temp = temp.split(' ')
            if temp[0:2] == ["MC","step"]:
                MCstep = int(temp[-1])
                #print(MCstep)
            if temp[0:2] == ["At", "point"]:
                #print(temp)
                site = temp[2].replace(")", "")
                site = site.replace("(", "")
                site = site.split(',')
                i = int(site[0])
                j = int(site[1])
                #print("({},{})".format(i,j))
                site_dict = {"step": MCstep, "i": i, "j": j, "line": n+start}
                df = pd.DataFrame([site_dict])
                catalog = pd.concat([catalog,df])
        self.site_catalog = catalog.reset_index()
        
    def check_for_missing_sites(self):
        len_file, testing_logfile = self._get_test_logfile()
        L = int(testing_logfile[5].split(' ')[-1])
        self.L = L
        nMC = int(testing_logfile[9].split(' ')[-1])
        self.nMC = nMC
        freq = int(testing_logfile[10].split(' ')[-1])
        self.freq = freq
        MC_steps = np.arange(0, nMC, freq)
        sites = self.site_catalog.copy()
        missing = []
        for n in MC_steps:
            sites_n = sites[sites["step"] == n]
            for i in range(0,L):
                sites_ni = sites_n[sites_n["i"] == i]
                for j in range(0,L):
                    sites_nij = sites_ni[sites_ni["j"] == j]
                    if sites_nij.empty:
                        missing.append([n,i,j])
                        site_dict = {"step": n, "i": i, "j": j, "line": np.nan}
                        df = pd.DataFrame([site_dict])
                        sites = pd.concat([sites, df])
        self.site_catalog = sites
        return missing
                    
    def _get_test_logfile(self):
        dir_path = None
        for item in os.listdir(self.test_run_dir):
            if item.startswith(self.dirheader):
                params=self.get_file_params(item,fordir=True) #if directory, not file, don't need to take off .csv
                if params==self.test_params: #find the right file
                    dir_path = self.test_run_dir+item #create path to that sub-directory
        if dir_path == None:
            self._error_message("No valid directory found in path {}".format(self.test_run_dir), "Params used:", self.test_params)
            sys.exit(1)
        else:        
            for file in os.listdir(dir_path):
                if file == self.test_log:
                    file_path = dir_path+"/"+file
                    f = open(file_path, "r")
                    testing_log = f.readlines()
                    len_file = int(len(testing_log))
                    f.close()
        return len_file, testing_log
        
                    
    def _search(self, file, phrase, testing = False):
        for n, line in enumerate(file):
            if phrase in line:
                if testing:
                    msg = "Found phrase at line "+str(n)
                    self._message(msg, phrase, line)
                return n
        self._message("Phrase not found in range.")
        return None
                
        
        
    
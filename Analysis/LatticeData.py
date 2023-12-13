#Authors: Casey Berger and Andy Esseln
#Last updated: 2023-12-13 by Casey 
#Last edit: added Andy's new functions, added _error_message and _message functions to clean up use of "print"

'''
Some proposed changes:

- Make multiple classes -- one for data utils and one for the data itself? Or maybe you need one for the correlation function data and one for the other data? The observables?
'''


import os
import shutil
import numpy as np
import pandas as pd



class LatticeData:
    def __init__(self, datadir = "/data/", header = "nonlinearsigma_data",
                 dirheader = "nlsigma_data", Gheader = "Gij_avg_nonlinearsigma_data", 
                 tol = 0.00001, palette = "viridis"):
        self.path = os.getcwd()+datadir #select location of data
        self.header = header #set the start of the filename for the data files
        self.dirheader = dirheader #set the start of the data directory name from the runs
        self.Gheader = Gheader #set the start of the filename for correlation function files
        self.tol = tol #set the error range for parameters -- this is for filtering
        self.palette = palette #option to change seaborn palette
        self.observables = ['Q_L', 'A_L', 'S_L', 'Xi_L'] #observables whose expectation values can be computed
        self.parameters = ["itheta", "beta", "length","nMC", "ntherm", "freq"] #parameters read in by the simulation code
        self.df_stats = pd.DataFrame() #initialize an empty dataframe to fill later
        self.df_stacked = False #set internal state for that dataframe
        
    #external functions / public
    def copy_data_from_directory(self, src_dir, dst_path = None):
        src_path = os.getcwd()+'/'+src_dir+'/'
        if dst_path is None:
            dst_path = self.path
        else:
            dst_path = dst_path
        for item in os.listdir(src_path):
            if item.startswith(self.dirheader):
                dir_path = src_path+item
                nMC = int(item.split("_")[-5])
                freq = int(item.split("_")[-1])
                for file in os.listdir(dir_path):
                    file_path = dir_path+"/"+file
                    if file.startswith(self.header):
                        f = open(file_path, "r")
                        len_file = int(len(f.readlines()))
                        len_complete = int(nMC/freq +1)
                        if len_file == len_complete:
                            shutil.copyfile(file_path, dst_path+file)
                        else:
                            status_msg = "Status: "+str(len_file)+" lines written"
                            self._message("run not yet complete",file[20:-4],status_msg)
                        f.close()
                    elif file.startswith(self.Gheader):
                        shutil.copyfile(file_path, dst_path+file)
                        
    def all_params(self):
        param_df = pd.DataFrame()
        files = self.get_data_files()
        for file in files:
            pdict = self.get_file_params(file)
            param_df = param_df.append(pdict,ignore_index=True)
        param_df["itheta/pi"] = param_df["itheta"]/np.pi
        return param_df
    
    def get_data(self,single_run = False, corr = False, suppress_output = True, reset_index = True,**kwargs):
        if single_run:
            missing_params = []
            for p in self.parameters:
                if p not in kwargs.keys():
                    missing_params.append(p)
            if len(missing_params) > 0:
                self._error_message("Missing parameters in input: ", missing_params)
                return None
        else:
            #self._message("Returning multiple runs")
            pass
        files = self.get_data_files(corr = corr)
        df = pd.DataFrame()
        for file in files:
            temp = pd.read_csv(file, skipinitialspace = True)
            pdict = self.get_file_params(file)
            if self.in_list(pdict, **kwargs):
                for key, value in pdict.items():
                    temp[key] = float(value)
                if not corr:
                    for observable in self.observables:
                        temp[observable+"_ta"]= self.ta(temp[observable])
                    F_py = self.calc_F(**pdict)
                    corr_length = self.calc_corr_length(temp["Xi_L"],temp["length"],F_py)
                    mass_gap = 1./corr_length
                    temp["corr_length_Re"] = corr_length.real
                    temp["corr_length_Im"] = corr_length.imag
                    temp["F_Re_py"] = F_py.real
                    temp["F_Im_py"] = F_py.imag
                    temp["mass_gap_Re"] = mass_gap.real
                    temp["mass_gap_Im"] = mass_gap.imag
                df = pd.concat([df,temp])
                if not suppress_output:
                    for key, value in pdict.items():
                        print(key, value)
        if reset_index:
            df.reset_index()
        return df
    
    def do_stats(self, therm = 0., stack = False, **kwargs):
        df = self.get_data(**kwargs)
        therm_condition = df["step"].astype(float) >= therm*df["nMC"].astype(float)
        df = df[therm_condition]
        df.drop(columns = ["step"], inplace = True)
        if stack:
            df_max = df.groupby(["length","itheta","beta","nMC","ntherm", "freq"]).max()
            df.drop(columns = ["dt"], inplace = True)
            df_means = df.groupby(["length","itheta","beta","nMC","ntherm", "freq"]).mean()
            ta_cols = [i+"_ta" for i in self.observables]
            df.drop(columns = ta_cols, inplace = True)
            df_sdevs = df.groupby(["length","itheta","beta","nMC","ntherm", "freq"]).std()
            df_all = df_means.join(df_sdevs,lsuffix = "_mean",rsuffix = "_std")
        else:
            df_max = df.groupby(["length","itheta","beta","nMC","ntherm","freq"]).max().reset_index()
            df.drop(columns = ["dt"], inplace = True)
            df_means = df.groupby(["length","itheta","beta","nMC","ntherm","freq"]).mean().reset_index()
            ta_cols = [i+"_ta" for i in self.observables]
            df.drop(columns = ta_cols, inplace = True)
            df_sdevs = df.groupby(["length","itheta","beta","nMC","ntherm","freq"]).std().reset_index()
            df_all = df_means.join(df_sdevs,lsuffix = "_mean",rsuffix = "_std")
            df_all.drop(columns = ["length_std", "itheta_std","beta_std", "nMC_std", "ntherm_std","freq_std"],inplace = True)
            df_all.rename(columns = {"length_mean":"length", "itheta_mean":"itheta","beta_mean":"beta", 
                           "nMC_mean":"nMC", "ntherm_mean":"ntherm", "freq_mean":"freq"}, inplace = True)
        df_all["time (sec)"] = df_max["dt"]
        df_all["time (min)"] = df_all["time (sec)"]/60.
        df_all["time (hr)"] = df_all["time (sec)"]/3600.
        self.df_stats = df_all
        self.df_stacked = stack
        return df_all
    
    def get_plot_data(self, obs = "Q_L", L = 10, beta = 1.6, nMC = 10000, ntherm = 1000, freq = 100):
        if len(self.df_stats) == 0:
            self._message("Generating dataframe of data with default statistical analysis")
            df = self.do_stats()
        df = self.df_stats
        if self.df_stacked:
            df.columns = pd.MultiIndex.from_product([df.columns, ["data"]])
            len_mask = df.index.get_level_values('length') == L 
            beta_mask = df.index.get_level_values('beta') == beta 
            nMC_mask = df.index.get_level_values('nMC') == nMC
            ntherm_mask = df.index.get_level_values('ntherm') == ntherm
            freq_mask = df.index.get_level_values('freq') == freq
            df = df[len_mask & beta_mask & nMC_mask & ntherm_mask & freq_mask]
            df = df.unstack(level = [0,2,3,4,5])
            df.columns = df.columns.droplevel(level = [1,2,3,4,5,6])
            x = df.index.to_numpy()
        else:
            len_mask = df['length'] == L 
            beta_mask = df["beta"] == beta 
            nMC_mask = df["nMC"] == nMC
            ntherm_mask = df["ntherm"] == ntherm
            freq_mask = df["freq"] == freq
            df = df[len_mask & beta_mask & nMC_mask & ntherm_mask & freq_mask]
            x = df["itheta"]
        y = df[obs+"_mean"]
        err = df[obs+"_std"]
        return x, y, err/np.sqrt(len(err))
    
    def get_corr_func(self,suppress_output = False,**kwargs):
        df = self.get_data(single_run = True, corr = True, suppress_output = suppress_output, **kwargs)
        length = kwargs["length"]
        df["i,j"] = df["i"]+df["j"]
        G_avg = df["G_avg"].to_numpy()
        G_avg = G_avg.reshape((length,length))
        return G_avg

    '''
    #old version (pre-Nov 8)
    def get_exceptional_configurations(self,src_dir):
        src_path = os.getcwd()+'/'+src_dir+'/' 
        config_df = pd.DataFrame() 
        for item in os.listdir(src_path):
            if item.startswith(self.dirheader): 
                dir_path = src_path+item 
                pdict = dict() 
                for file in os.listdir(dir_path): 
                    file_path = dir_path+"/"+file
                    if file.startswith(self.header): 
                        pdict = self.get_file_params(file) 
                for file in os.listdir(dir_path): 
                    if file.startswith("config"): 
                        config_dict = dict() #create empty dictionary to store config number and number of exceptional sites in that config
                        config_dict.update(pdict) #add parameter dictionary to the config dict
                        config_num = int(file[0:-4].split("_")[-1]) #pull the number from the filename
                        config_dict["config"] = config_num #store the config number
                        temp = pd.read_csv(file_path, skipinitialspace = True) #read the config file
                        num_N = temp['exceptional'].value_counts()['N'] #count all the "N"s to avoid key error            
                        num_exc = len(temp['exceptional']) - num_N #number of exc
                        config_dict["num_exc"] = num_exc #add to dictionary
                        config_df = config_df.append(config_dict,ignore_index=True)
        config_df["any_exc"] = config_df["num_exc"]>0 #flag all configurations that have any exceptional sites
        config_df["any_exc"] = config_df["any_exc"].astype(int) #store as 0s and 1s instead of bool for counting
        return config_df
    '''

    def get_exceptional_configurations(self,src_dir):
        src_path = os.getcwd()+'/'+src_dir+'/' #create path to run directory
        config_df = pd.DataFrame() #create empty data frame
        for item in os.listdir(src_path):#loop over all files in the run directory
            if item.startswith(self.dirheader): #pick out the sub-directories (each is an individual run of the code)
                dir_path = src_path+item #create path to the sub-directory
                pdict = dict() #create empty dictionary for parameters
                nMC = int(item.split("_")[-5])
                freq = int(item.split("_")[-1])
                for file in os.listdir(dir_path): #loop over every file in the run sub-directory
                    file_path = dir_path+"/"+file #create path to the file we're looking at
                    if file.startswith(self.header): #if it's the full observable logfile
                        f = open(file_path, "r")
                        len_file = int(len(f.readlines()))
                        len_complete = int(nMC/freq +1)
                        if len_file != len_complete:
                            status_msg = "Status: "+str(len_file)+" lines written"
                            self._message("run not yet complete",file[20:-4],status_msg)
                            good=False
                        else:
                            pdict = self.get_file_params(file) #add the parameters for this run to pdict
                            good=True    
                        f.close()
                if good==True:
                    for file in os.listdir(dir_path): #loop again over every file
                        if file.startswith("config"): #find just the config files
                            file_path = dir_path+"/"+file #create path to the file we're looking at
                            config_dict = dict() #create empty dictionary to store config number and number of exceptional sites in that config
                            config_dict.update(pdict) #add parameter dictionary to the config dict
                            config_num = int(file[0:-4].split("_")[-1]) #pull the number from the filename
                            config_dict["config"] = config_num #store the config number
                            temp = pd.read_csv(file_path, skipinitialspace = True) #read the config file
                            num_N = temp['exceptional'].value_counts()['N'] #count all the "N"s to avoid key error            
                            num_exc = len(temp['exceptional']) - num_N #number of exc
                            config_dict["num_exc"] = num_exc #add to dictionary
                            config_df = config_df.append(config_dict,ignore_index=True)
        config_df["any_exc"] = config_df["num_exc"]>0 #flag all configurations that have any exceptional sites
        config_df["any_exc"] = config_df["any_exc"].astype(int) #store as 0s and 1s instead of bool for counting
        return config_df
    
    def find_exc(self,src_dir,**kwargs): 
        '''
        takes a dictionary of parameters as kwargs just like get_data for one run
        src_dir is the raw data directory, not the processed directory

        #note from Andy:
        #to make graphs with this function, I do something like...
        corr_params = {"itheta": 0.5, "beta": 1.6,"length": 40,"nMC": 50000, "ntherm": 0, "freq": 100}
        excarray=analyzer.find_exc('/data/run_2023_10_26_exc_config_test',**corr_params)
        aconfig=excarray[:,:,0] #0 if you want to show the first configuration, or put in another number
        plt.imshow(aconfig)
        '''
        missing_params=[]
        for p in self.parameters:
            if p not in kwargs.keys():
                missing_params.append(p)
        if len(missing_params) > 0:
            self._error_message("Missing parameters in input: ", missing_params)
            return None #tells you if you didn't put in enough paramters, again just like get_data
        src_path = os.getcwd()+'/'+src_dir+'/' #create path to run directory
        for item in os.listdir(src_path):#loop over all files in the run directory
            if item.startswith(self.dirheader): #pick out the sub-directories (each is an individual run of the code)
                params=self.get_file_params(item,fordir=True) #if directory, not file, don't need to take off .csv
                if params==kwargs: #find the right file
                    dir_path = src_path+item #create path to that sub-directory
                    length=params['length']
                    freq=params['freq']
                    nMC=params['nMC']
                    for file in os.listdir(dir_path): #loop over every file in the run sub-directory
                        file_path = dir_path+"/"+file #create path to the file we're looking at
                        if file.startswith(self.header): #if it's the full observable logfile
                            f = open(file_path, "r")
                            len_file = int(len(f.readlines()))-1 #take out top line
                            len_complete = int(nMC/freq)
                            if len_file != len_complete: #we can still look at where the exceptional configurations are in an imcomplete file!
                                status_msg = "Status: "+str(len_file)+" lines written. Continuing..."
                                self._message("run not yet complete",file[20:-4],status_msg)
                    excarray=np.empty((length,length,len_file),dtype=bool) #3D array - rows in lattice, columns in lattice, number of (completed) configs
                    for file in os.listdir(dir_path): #loop over every file in the subdirectory
                        file_path = dir_path+"/"+file #make path to that file
                        if file.startswith("config"): #pick out just the config files
                            cfn=int(file[7:-4]) #get number of the configuration from file name
                            temp=pd.read_csv(file_path, skipinitialspace = True)
                            for row in range(length):
                                yns=temp['exceptional'][length*row:length*(row+1)] #get data for first row
                                excarray[row,:,int(cfn/freq)]=self._yn2bool(yns) #turn it into booleans, put it in the array...
                                #1st index will count up as we go through this loop, go through each row
                                #2nd index is : because you're filling in the whole row
                                #3rd index puts it in the right place for the configuration number - ie, for config330 with frequency=10, it goes in the 33rd
        return excarray

    
    def convert_config_spherical(self,src_dir):
        '''
        Add Andy's to_spherical function in here to convert each config to spherical??
        '''
        count = 0 #testing
        src_path = os.getcwd()+'/'+src_dir+'/' 
        config_df = pd.DataFrame() 
        for item in os.listdir(src_path):
            if item.startswith(self.dirheader): 
                dir_path = src_path+item 
                pdict = dict() 
                for file in os.listdir(dir_path): 
                    file_path = dir_path+"/"+file
                    if file.startswith(self.header): 
                        pdict = self.get_file_params(file) 
                for file in os.listdir(dir_path): 
                    if file.startswith("config"): 
                        df = pd.read_csv(file)
                        count +=1
                    if count > 0:
                        break
    
    def to_spherical(x,y,z):
        '''
        note from Andy:
        if z=1 - phi points straight up - then the polar angle returns undefined. That makes sense, but 
        is it better to assign it a nan in this situation rather than getting it via the divide by 0 error? 
        Or to assign a set value (like 0?) to the polar angle in that case?

        It was such a big deal with the other code that I tried to avoid using arccos's... probably for no 
        reason, but it works this way too so might as well

        I thought about putting in a line to check that the components you put have length one and to give 
        an error if they don't, but I figured that would be checked elsewhere in the code, and also that 
        rounding somewhere might cause them to add to almost-but-not-quite 1, and so including a line like 
        that would be making problems where none exist
        '''
        azimuthal=np.arcsin(np.sqrt(1-z**2))
        polar=np.arcsin(y/(np.sqrt(1-z**2)))
        if x<0:
            polar=np.pi-polar
        return polar,azimuthal
        
    #internal functions / private    
    
    def _yn2bool(self,series):
        #required for find_exc()
        boollist=[]
        for i in series:
            if i=='Y':
                boollist.append(True)
            elif i=='N':
                boollist.append(False)
            else:
                self._message('not Y/N')
        return boollist
    
    def get_data_files(self, corr = False):
        data_files = []
        file_header = self.header
        if corr:
            file_header = self.Gheader
        for file in os.listdir(self.path):
            if file.startswith(file_header):
                data_files.append(self.path+file)
        return data_files

    def get_file_params(self, file, fordir = False):
        if not fordir:
            file = file[:-4]#remove ".csv" before splitting
        temp = file.split("_")
        pdict = dict()
        pdict[temp[-2]] = int(temp[-1]) #frequency
        pdict[temp[-4]] = int(temp[-3]) #nMc
        pdict[temp[-6]] = int(temp[-5]) #ntherm
        pdict[temp[-8]] = float(temp[-7]) #itheta
        pdict[temp[-10]] = float(temp[-9]) #beta
        pdict["length"] = int(temp[-11]) #length
        return pdict
    
    def in_list(self,pdict,**kwargs):
        count = 0
        for key, value in kwargs.items():
            if (pdict[key] > value+self.tol) or (pdict[key] <  value-self.tol):
                break
            count += 1
        if count == len(kwargs.keys()):
            return True
        else:
            return False
        
    def calc_F(self, **kwargs):
        #fix this to actually get the 2-point correlation function instead of the average G
        G_avg = self.get_corr_func(suppress_output = True, **kwargs)
        L = kwargs["length"]
        F = complex(0,0)
        i = complex(0,1)
        for x1 in range(0,L):
            for x2 in range(0,L):
                F += (np.exp(2.*np.pi*i*x1/L)+ np.exp(2.*np.pi*i*x2/L))*G_avg[x1,x2]
        return 0.5*F
    
    def calc_corr_length(self,Xi,L,F_py):
        Xi = Xi.to_numpy()
        L = L.to_numpy()
        return np.sqrt(Xi/F_py)/(2.*np.sin(np.pi/L))
    
    def next_pow_two(self,n):
        i = 1
        while i < n:
            i = i << 1
        return i

    def autocorr_func_1d(self, x, norm=True):
        '''
        This function comes straight from https://dfm.io/posts/autocorr/ and computes the 
        autocorrelation of a MCMC
        
        You feed into it x, which is the Markov Chain
        '''
        x = np.atleast_1d(x)
        if len(x.shape) != 1:
            raise ValueError("invalid dimensions for 1D autocorrelation function")
            
        n = self.next_pow_two(len(x))
        
        # Compute the FFT and then (from that) the auto-correlation function
        f = np.fft.fft(x - np.mean(x), n=2*n)
        acf = np.fft.ifft(f * np.conjugate(f))[:len(x)].real
        acf /= 4*n
        
        # Optionally normalize
        if norm:
            acf /= acf[0]
        return acf
    
    def autocorrelation(self, df, norm = True):
        for observable in self.observables:
            df["acf_"+observable] = self.autocorr_func_1d(df[observable], norm=norm)
        return df
    
    def ta(self,data_array):
        acf = self.autocorr_func_1d(data_array, norm=False)#switch back to true
        decorr = np.where(acf < 0.3)
        return decorr[0][1]
    
    def _error_message(self, *args):
        print("Error: ")
        for arg in args:
            print(arg)
            
    def _message(self, *args):
        for arg in args:
            print(arg)
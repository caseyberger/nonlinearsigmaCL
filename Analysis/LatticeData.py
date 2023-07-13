import os
import numpy as np
import pandas as pd

class LatticeData:
    def __init__(self, datadir = "/data/", header = "nonlinearsigma_data",
                 Gheader = "Gij_avg_nonlinearsigma_data", tol = 0.00001,
                 palette = "viridis"):
        self.path = os.getcwd()+datadir #select location of data
        self.header = header #set the start of the filename for the data files
        self.Ghead = Gheader
        self.tol = tol
        self.palette = palette #option to change seaborn palette
        self.observables = ['Q_L', 'A_L', 'S_L', 'Xi_L']
        self.parameters = ["itheta", "beta", "length","nMC", "ntherm"]
        
    def get_data_files(self, corr = False):
        data_files = []
        file_header = self.header
        if corr:
            file_header = self.Ghead
        for file in os.listdir(self.path):
            if file.startswith(file_header):
                data_files.append(self.path+file)
        return data_files

    def get_file_params(self, file):
        file = file[:-4]
        temp = file.split("_")
        pdict = dict()
        pdict[temp[-2]] = int(temp[-1]) #nMc
        pdict[temp[-4]] = int(temp[-3]) #ntherm
        pdict[temp[-6]] = float(temp[-5]) #itheta
        pdict[temp[-8]] = float(temp[-7]) #beta
        pdict["length"] = int(temp[-9]) #length
        return pdict
    
    def all_params(self):
        param_df = pd.DataFrame()
        files = self.get_data_files()
        for file in files:
            pdict = self.get_file_params(file)
            param_df = param_df.append(pdict,ignore_index=True)
        return param_df
    
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
    
    def get_data(self,single_run = False, corr = False, suppress_output = True, **kwargs):
        if single_run:
            missing_params = []
            for p in self.parameters:
                if p not in kwargs.keys():
                    missing_params.append(p)
            if len(missing_params) > 0:
                print("Missing parameters in input: ")
                print(missing_params)
                return None
        else:
            #print("Returning multiple runs")
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
        return df
    
    def do_stats(self, therm = 0., **kwargs):
        df = self.get_data(**kwargs)
        therm_condition = df["step"].astype(float) >= therm*df["nMC"].astype(float)
        df = df[therm_condition]
        df.drop(columns = ["step"], inplace = True)
        df_max = df.groupby(["length","itheta","beta","nMC","ntherm"]).max()
        df.drop(columns = ["dt"], inplace = True)
        df_means = df.groupby(["length","itheta","beta","nMC","ntherm"]).mean()
        ta_cols = [i+"_ta" for i in self.observables]
        df.drop(columns = ta_cols, inplace = True)
        df_sdevs = df.groupby(["length","itheta","beta","nMC","ntherm"]).std()
        df_all = df_means.join(df_sdevs,lsuffix = "_mean",rsuffix = "_std")
        df_all["time (sec)"] = df_max["dt"]
        df_all["time (min)"] = df_all["time (sec)"]/60.
        df_all["time (hr)"] = df_all["time (sec)"]/3600.
        return df_all
    
    def get_corr_func(self,suppress_output = False,**kwargs):
        df = self.get_data(single_run = True, corr = True, suppress_output = suppress_output, **kwargs)
        length = kwargs["length"]
        df["i,j"] = df["i"]+df["j"]
        G_avg = df["G_avg"].to_numpy()
        G_avg = G_avg.reshape((length,length))
        return G_avg
    
    def calc_F(self, **kwargs):
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
            df["acf_"+observable] = self.autocorr_func_1d(df[observable], norm)
        return df
    
    def ta(self,data_array):
        acf = self.autocorr_func_1d(data_array, norm=True)
        decorr = np.where(acf < 0.3)
        return decorr[0][1]
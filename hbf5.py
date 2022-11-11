import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess

def bin_log_lin(mn, mx, nb_log, nb_lin):
    """ Takes an interval and splits it into layered bins.
        First layer is split logarithmically, second is linearly.
        Arguments:
            mn = lower bound of interval
            mx = upper bound of interval
            nb_log = number of logarithmic bins
            nb_lin = numer of linear bins
        Returns:
            bins = returns list of bin intervals
    """
    
    # create logarithmic bins
    log_bins = np.logspace(np.log10(mn), np.log10(mx), nb_log+1)

    print(log_bins)

    # create linear bins
    bins = np.linspace(log_bins[:-1], log_bins[1:], nb_lin)

    print(bins)

    # transpose and flatten
    bins = bins.transpose().flatten()

    return bins

#-----------------------------------------------------------------------------------------

FB_PATH = "/home/tomas/Documents/cmu/research/fewbody/fewbody-pn/" # directory fewbody is in
FB_PROCESS = "binsingle_hbf5"    # fewbody process to be used
DATA_PATH = '/home/tomas/Documents/cmu/research/hbf5/data/'
FIGS_PATH = '/home/tomas/Documents/cmu/research/hbf5/figs/'
N = 500                 # number of interactions to be computed
a = 10.0                # semi-major axis of initial binary (in AU)
ecc = 0.0               # eccentricity of initial binary
hb_C = 4.0              # C parameter from Hut & Bahcall 1983
hb_D = 0.6*(1+ecc)      # D parameter from Hut & Bahcall 1983
v_min = 0.5             # minimum of velocity sampling range (in units of critical velocity)
v_max = 16.0            # maximum of velocity sampling range ""
n_vbins_log = 5         # number of logarithmic velocity bins
n_vbins_lin = 20        # number of linear velocity bins per log bin

# make parameter bins
v_bounds = bin_log_lin(v_min, v_max, n_vbins_log, n_vbins_lin)
print(v_bounds)

# make any necessary files
subprocess.run("make -C "+FB_PATH, shell=True)

# generate array to store initial conditions and classification
dat_list = []
weird_no = []

# iterate and integrate
for i in range(N):

    # generate three random numbers: one for velocity, one for impact parameter, and one for fewbody seed
    rn_set = np.random.rand(3)

    # generate velocity, iterating over bins
    bin_no = i%(n_vbins_log*n_vbins_lin)
    v = v_bounds[bin_no]
    #v = v_bounds[bin_no] + rn_set[0]*(v_bounds[bin_no+1]-v_bounds[bin_no]) 

    # generate impact parameter using velocity (in units of semi-major axis)
    b_max = (hb_C/v+hb_D)
    b = b_max*(rn_set[1])**0.5

    # calculate cross-section (in units of orbital area)
    #sigma_frac = (b_max)**2

    # generate seed
    s = int(rn_set[2]*2**32)
    
    # run process; "2> /dev/null" directs (>) stderr output (2) to the trash vent of linux (/dev/null)
    binout = subprocess.check_output(FB_PATH+FB_PROCESS+" --vinf "+str(v)+" --b "+str(b)
                                     +" --seed "+str(s)+" 2> /dev/null", shell=True)
    mess = binout.decode("utf-8")
    clss = mess[-2]
    if clss == 'i':
        dat_list.append([v,b,s,0])
    elif clss == 'e':
        dat_list.append([v,b,s,1])
    elif clss == 'r': 
        dat_list.append([v,b,s,2])
    elif clss == 'u':
        print(i," WARNING: type is not classified!")
        weird_no.append(i)

    if (i/N*100)%1==0:
        print('%.2f percent done (i=%d)' % (i/N*100, i))

# wrap output into pandas dataframe
dat = pd.DataFrame(dat_list, columns=['v','b','s','clss'])
print(dat)

# generate count matrix (velocity x class)
counts = pd.DataFrame(index=dat.v.unique())
counts['ionization'] = dat[dat.clss==0].groupby('v')['b'].count()
counts['exchange'] = dat[dat.clss==1].groupby('v')['b'].count()
counts['resonance'] = dat[dat.clss==2].groupby('v')['b'].count()
print(counts)

# generate cross section matrix
xcs = counts.mul((hb_C/counts.index+hb_D)**2, axis=0)
xcs.to_csv(DATA_PATH+"xcs_%s_N-%d_e-%.3f.csv" % (FB_PROCESS, N, ecc))
print(xcs)

# check for non-classified interactions
if len(weird_no) != 0:
    print("WARNING: ",len(weird_no)," non-classifiable interactions were found ",weird_no)

for i in xcs.columns:
    plt.scatter(xcs.index,xcs[i],label=i)
plt.title('N=%d, e=%f' % (N,ecc))
plt.xlabel('velocity [v_crit]')
plt.ylabel('cross section [orbital area]')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig("figs/%s_N-%d_e-%.3f.png" % (FB_PROCESS, N, ecc))
plt.close()

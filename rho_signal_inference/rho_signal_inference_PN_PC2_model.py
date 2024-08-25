# -*- coding: utf-8 -*-
"""
Rho inference
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
from sklearn.decomposition import PCA
import seaborn as sns

##############################################################################
############### Utility functions
##############################################################################

def set_font_sizes(SMALL_SIZE=14, MEDIUM_SIZE=16, LARGE_SIZE=20):
    '''
    Sets font size for matplotlib
    From: https://stackoverflow.com/a/39566040
    '''
    font = {#'family':'sans-serif',
            'sans-serif':['Arial'],
            'size': SMALL_SIZE}
    plt.rc('font', **font)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=LARGE_SIZE)   # fontsize of the figure title

set_font_sizes()

np.random.seed(1234)

# set up ORN and PN colors
d_color = {'ORN': '#ffc951',
           'PN':  '#8621f1'}

def generate_correlated_x(x, rho):
    '''
    Given a sample of numbers x, 
    and a correlation rho,
    returns sample of numbers with correlation rho to x
    '''
    z = np.random.normal(size=len(x))
    return rho * x + np.sqrt(1 - rho**2) * z

def get_r(x, y):
    return np.corrcoef(x, y)[0, 1]


def run_forward_analysis(rho, rpred, rb, Nflies):
    '''
    Given a biological signal strength rho, 
    a pearson correlation coefficient for calcium / brp repeatability rpred,
    a pearson correlation coefficient for behavior repeatability rb,
    and a number of flies,
    runs a simulation propagating from latent calcium / brp to a 
    final R^2 between observed calcium /brp and behavior (R^2_c/brp,b)
    '''
    
    # simulate latent calcium / brp
    x_pred_latent = np.random.normal(size=Nflies)
    # generate behavior through rho
    x_behav_latent = generate_correlated_x(x_pred_latent, rho)
    
    # simulate observed calcium / brp
    x_pred_obs1 = generate_correlated_x(x_pred_latent, rpred)
    x_pred_obs2 = generate_correlated_x(x_pred_latent, rpred)

    # simulate observed behavior
    x_behav_obs1 = generate_correlated_x(x_behav_latent, rb)
    x_behav_obs2 = generate_correlated_x(x_behav_latent, rb)
    
    r_squared_predb = np.corrcoef(x_pred_obs1, x_behav_obs1)[0,1]**2
    
    return r_squared_predb, x_pred_latent, x_behav_latent, x_pred_obs1, x_pred_obs2, x_behav_obs1, x_behav_obs2
    

##############################################################################
############### Data
##############################################################################


df_persistence = pd.read_csv('data/OCT-MCH_behavior_persistence.csv')
############## BEHAVIOR - BEHAVIOR
OCTMCH_h0 = df_persistence['OCT-MCH hour 0'].values 
OCTMCH_h3 = df_persistence['OCT-MCH hour 3'].values

n_behav = len(OCTMCH_h0)


# get R^2
r_bb = get_r(OCTMCH_h0, OCTMCH_h3)
R2_bb = r_bb**2

# perform bootstraps
nboot = 10000
R2_bb_bootstraps = np.zeros(nboot)
for i in range(nboot):
    bootstrap_indices = np.random.choice(np.arange(n_behav), size=n_behav, replace=True)
    x_bb_bootstrap = OCTMCH_h0[bootstrap_indices]
    y_bb_bootstrap = OCTMCH_h3[bootstrap_indices]
    R2_bb_bootstraps[i] = get_r(x_bb_bootstrap, y_bb_bootstrap)**2


#
# plot results
#
fig, axs = plt.subplots(1, 2, figsize=(10,6))

# original data
axs[0].set_title('OCT-MCH preference persistence (3 h)'+'\n'+\
                 r'$n=$ ' + f'{n_behav}, ' + r'$R^2=$' + '{:.2f}'.format(R2_bb))
axs[0].scatter(OCTMCH_h0, OCTMCH_h3, color='k')
axs[0].set_xlabel('odor preference (t = 0 h)')
axs[0].set_ylabel('odor preference (t = 3 h)')
axs[0].set_xlim(-0.8, 0.5)
axs[0].set_ylim(-0.8, 0.5)

# result from bootstrapping
axs[1].set_title(f'{nboot} bootstrap replicates')
axs[1].hist(R2_bb_bootstraps, bins=50, color='0.6',)
axs[1].axvline(R2_bb, c='k', ls='--', label=r'$R^2=$' + '{:.2f}'.format(R2_bb))
axs[1].set_ylabel('number of bootstraps')
axs[1].set_xlabel(r'$R^2$ of bootstrapped data')
axs[1].legend()

plt.suptitle(r'bootstrap analysis of behavior-behavior $R^2$'+'\n(figure 1 - figure supplement 1 D)')
plt.tight_layout()
plt.show()

############## CALCIUM - BEHAVIOR

df_PN_PCs = pd.read_csv('data/PN_PC2_interp_scores.csv')

PN_PC2_scores = df_PN_PCs['PN PC2 scores (Fig 1M)'].values
PN_OCT_MCH_behav_scores = df_PN_PCs['OCT-MCH behavior scores'].values 


PNPC2_OCTMCH_pred = PN_PC2_scores
PNPC2_OCTMCH_meas = PN_OCT_MCH_behav_scores

n_cb = len(PNPC2_OCTMCH_pred)

# get R^2
r_cb = get_r(PNPC2_OCTMCH_pred, PNPC2_OCTMCH_meas)
R2_cb = r_cb**2

# perform bootstraps
nboot = 10000
R2_cb_bootstraps = np.zeros(nboot)
for i in range(nboot):
    bootstrap_indices = np.random.choice(np.arange(n_cb), size=n_cb, replace=True)
    x_cb_bootstrap = PNPC2_OCTMCH_pred[bootstrap_indices]
    y_cb_bootstrap = PNPC2_OCTMCH_meas[bootstrap_indices]
    R2_cb_bootstraps[i] = get_r(x_cb_bootstrap, y_cb_bootstrap)**2


#
# plot results
#
fig, axs = plt.subplots(1, 2, figsize=(10,6))

# original data
axs[0].set_title('OCT-MCH preference predicted by PN Ca++ PC2'+'\n'+\
                 r'$n=$ ' + f'{n_cb}, ' + r'$R^2=$' + '{:.2f}'.format(R2_cb))
axs[0].scatter(PNPC2_OCTMCH_pred, PNPC2_OCTMCH_meas, color=d_color['PN'])
axs[0].set_xlabel('predicted preference (z-scored)')
axs[0].set_ylabel('measured preference (z-scored)')
axs[0].set_xlim(-3, 2.3)
axs[0].set_ylim(-3, 2.3)

# result from bootstrapping
axs[1].set_title(f'{nboot} bootstrap replicates')
axs[1].hist(R2_cb_bootstraps, bins=np.arange(0, 0.551, 0.01), color='0.6',)
axs[1].axvline(R2_cb, c='k', ls='--', label=r'$R^2=$' + '{:.2f}'.format(R2_cb))
axs[1].set_ylabel('number of bootstraps')
axs[1].set_xlabel(r'$R^2$ of bootstrapped data')
axs[1].legend()

plt.suptitle(r'bootstrap analysis of calcium-behavior $R^2$'+'\n(figure 1 M)')
plt.tight_layout()
plt.show()



############## CALCIUM - CALCIUM

PNdata = pd.read_csv('data/gh146flyaverage.txt', header=None).values


nPN = PNdata.shape[1]
ndim = 65

PNpca = PCA(); PNpca.fit(PNdata.T)

PNshuffled = np.zeros(PNdata.shape)
for j in range(nPN):
    PNshuffled[:, j] = PNdata[np.random.permutation(ndim), j]
     
    
PNshuffpca = PCA()
PNshuffpca.fit(PNshuffled.T)



diffratPN = PNpca.explained_variance_ratio_ - PNshuffpca.explained_variance_ratio_
PNcrossi = np.where(diffratPN < 0)[0][0]
R2_cc = np.sum(PNpca.explained_variance_ratio_[:PNcrossi])


def perform_PN_boot():
    # perform bootstrapping
    bootstrapPNis = np.random.choice(np.arange(nPN), nPN, replace=True)    
    bootstrapPNdata = PNdata[:,bootstrapPNis]

    # shuffle the bootstrapped PN data
    bootstrapPNshuffled = np.zeros(bootstrapPNdata.shape)
    for j in range(nPN):
        bootstrapPNshuffled[:, j] = bootstrapPNdata[np.random.permutation(ndim), j]
        
    # perform PCA
    PNbootpca = PCA()
    PNbootpca.fit(bootstrapPNdata.T)
    PNbootshuffpca = PCA()
    PNbootshuffpca.fit(bootstrapPNshuffled.T)

    # find summed variance for PCs where var explained > that of shuffled
    diffratPN = PNbootpca.explained_variance_ratio_ - PNbootshuffpca.explained_variance_ratio_
    PNcrossi = np.where(diffratPN < 0)[0][0]
    PN_explained_var = np.sum(PNbootpca.explained_variance_ratio_[:PNcrossi])
    
    return PN_explained_var


nboot = 1000
all_PN_explained_vars = np.zeros(nboot)

for i in range(nboot):
    all_PN_explained_vars[i] = perform_PN_boot()




#
# plot results
#
ndimplot = 10

fig, axs = plt.subplots(1, 2, figsize=(10,6))

# original data
axs[0].set_title('variance explained by PN calcium\ncompared to shuffled control'+'\n'+\
                 r'$n=$ ' + f'{nPN}, ' + r'$R^2=$' + '{:.2f}'.format(R2_cc))
axs[0].plot(1+np.arange(ndimplot), PNpca.explained_variance_ratio_[:ndimplot], label='PN calcium')
axs[0].plot(1+np.arange(ndimplot), PNshuffpca.explained_variance_ratio_[:ndimplot], label='shuffled PN calcium')
axs[0].axvline(PNcrossi, c='k', label='cutoff PC')
axs[0].legend()
axs[0].set_xlabel('principal component')
axs[0].set_ylabel('fraction variance explained')
axs[0].set_xticks(1+np.arange(ndimplot))


# result from bootstrapping
axs[1].set_title(f'{nboot} bootstrap replicates')
axs[1].hist(all_PN_explained_vars, bins=50, color='0.6',)
axs[1].axvline(R2_cc, c='k', ls='--', label=r'$R^2=$' + '{:.2f}'.format(R2_cc))
axs[1].set_ylabel('number of bootstraps')
axs[1].set_xlabel(r'$R^2$ of bootstrapped data')
axs[1].legend()

plt.suptitle(r'bootstrap analysis of calcium-calcium $R^2$')
plt.tight_layout()
plt.show()





##############################################################################
############### Run simulations
##############################################################################


N_exps = 10000
my_rhos = np.linspace(0, 1, 101)

N_rhos = len(my_rhos)

calculated_rpredictorb2s = np.zeros((N_exps, N_rhos))

check_predictor_rsquareds = []
check_behav_rsquareds = []

for ii in range(N_exps):
    if (ii % 500 == 0) & (ii > 0):
        print(f'{ii} experiments done')
    for rr in range(N_rhos):
        
        ## bootstrap behavior - behavior
        behav_bootstrap_indices = np.random.choice(np.arange(n_behav), size=n_behav, replace=True)
        cur_r2b = get_r(OCTMCH_h0[behav_bootstrap_indices], 
                        OCTMCH_h3[behav_bootstrap_indices])**2
        r_b = cur_r2b**(1/4)
               
        ## bootstrap calcium - calcium
        r_predictor = perform_PN_boot()**(1/4)
        
        
        cur_rho = my_rhos[rr]
        
        # run forward results
        r_squared_predictorb, x_predictor_latent, x_behav_latent, \
            x_predictor_obs1, x_predictor_obs2, x_behav_obs1, x_behav_obs2 \
                = run_forward_analysis(cur_rho, r_predictor, r_b, n_cb)
    
        check_predictor_rsquareds.append(np.corrcoef(x_predictor_obs1, x_predictor_obs2)[0,1]**2)
        check_behav_rsquareds.append(np.corrcoef(x_behav_obs1, x_behav_obs2)[0,1]**2)
        
        calculated_rpredictorb2s[ii, rr] = r_squared_predictorb
        
#############


predictor_color = d_color['PN']




################################

# check that R^2 between observed calcium/behavior matches empirical R^2
fig, axs = plt.subplots(1, 2, figsize=(10,5), sharey=True)
axs[0].hist(check_predictor_rsquareds, color=predictor_color, bins=50, alpha=0.5, density=True)
axs[0].axvline(r_predictor**4, c='k', ls='--')
#axs[0].set_title('{} repeatability check'.format(predictor_type))
axs[0].set_xlabel('predictor $R^2$')
axs[0].set_ylabel('# simulations')
axs[1].hist(check_behav_rsquareds, color=predictor_color, bins=50, alpha=0.5, label='simulations', density=True)
axs[1].axvline(r_b**4, c='k', ls='--', label='empirical')
axs[1].set_title('behavior repeatability check')
axs[1].set_xlabel('behavior $R^2$')
axs[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
#axs[1].hist(r2s_bootstrap_manual**(1/4), bins=50, density=True, color='k', alpha=0.5)
plt.show()



# count up simulations with given rhos

cb_hist_bins = np.arange(0, 0.6001, 0.01)
cb_hist_cnts, fo_ = np.histogram(R2_cb_bootstraps, bins=cb_hist_bins)
cb_hist_fracs = cb_hist_cnts / np.sum(cb_hist_cnts)
cb_hist_centers = cb_hist_bins[:-1] + (cb_hist_bins[1:] - cb_hist_bins[:-1])/2

n_hist_bins = len(cb_hist_cnts)


cnt_arr = np.zeros((n_hist_bins, len(my_rhos)))

for i in range(n_hist_bins):
    cur_lb = cb_hist_bins[i]
    cur_ub = cb_hist_bins[i+1]

    sim_cnt_per_rho = ((calculated_rpredictorb2s >= cur_lb) & (calculated_rpredictorb2s < cur_ub)).sum(0)
    
    cnt_arr[i, :] = sim_cnt_per_rho
    

resvals = (cnt_arr.T * cb_hist_fracs).T.sum(0)
resvals /= np.trapz(resvals, x=my_rhos)

drho = my_rhos[1] - my_rhos[0]

my_cdf = np.cumsum(resvals)/sum(resvals)


# quantiles
qs = [0.05, 0.5, 0.95]
for q in qs:
    where_cutoff = np.where(my_cdf <= q)[0][-1]
    rho_at_cutoff = my_rhos[where_cutoff]
    print('{:0>2} percentile: rho = {:.2f}'.format(q, rho_at_cutoff))
    
    
### SUP FIG

fig, axs = plt.subplots(2, 2, figsize=(6, 6), 
                        width_ratios = [5, 1], 
                        height_ratios=[5,1], 
                        sharey='row',
                        sharex='col')#; axs = [axs]

cmap = colormaps['binary']

plot_maxx = 0

ax0 = axs[0,0]


cur_median_rpredictorb2s = np.median(calculated_rpredictorb2s, 0)
for i in range(1, 10):
    cur_pctile = i/10 - 1/20
    c_arg = 1 - np.abs(cur_pctile - 0.5 + 0.05)*2
    my_clr = cmap(c_arg)
    
    cur_pctile_rsquared_maxes = np.quantile(calculated_rpredictorb2s, cur_pctile, axis=0)
    next_pctile_rsquared_maxes = np.quantile(calculated_rpredictorb2s, cur_pctile + 0.1, axis=0)        
    
    if i < 10:
        ax0.fill_betweenx(my_rhos, 
                         x1 = cur_pctile_rsquared_maxes, 
                         x2 = next_pctile_rsquared_maxes, 
                         color=my_clr, lw=0, alpha=0.55)
    
        maxval = np.max(np.concatenate((cur_pctile_rsquared_maxes, next_pctile_rsquared_maxes)))
        #print(i, maxval)
        plot_maxx = np.max([plot_maxx, maxval])

ax0.plot(cur_median_rpredictorb2s, my_rhos, c='k')#, label='median', c='k')
#cur_ax.legend(loc='upper left', bbox_to_anchor=(1.03, 1), borderaxespad=0)

ax0.set_ylim(0, 1)

for ax in [ax0]:
    ax.set_xlim(0, plot_maxx)


dspace = 0.2
#plt.subplots_adjust(hspace=0.3, wspace=0.1)
plt.subplots_adjust(hspace=dspace, wspace=dspace) 

ax1 = axs[1, 0]

sns.kdeplot(R2_cb_bootstraps, ax=ax1, color=predictor_color, lw=2,
            alpha=0.4, fill=True)
ax1.axvline(R2_cb, c='k', ls='--', label=r'$R^2=$' + '{:.2f}'.format(R2_cb))

ax0.set_ylabel(r'$\rho_{signal}$')
ax1.set_xlabel(r'$R^2_{c,b}$')
ax1.set_ylabel('')
ax1.set_yticks([])

ax2 = axs[0,1]

ax2.set_xticks([])

axs[1,1].axis('off')


ax2.plot(resvals, my_rhos, color=predictor_color, lw=2)
ax2.fill_betweenx(my_rhos, 0, resvals, color=predictor_color, alpha=0.4)
ax2.set_xlim(0, 1.05*max(resvals))

qs = [0.05, 0.5, 0.95]


for q in qs:
    where_cutoff = np.where(my_cdf <= q)[0][-1]
    rho_at_cutoff = my_rhos[where_cutoff]
    #print('{:0>2} percentile: rho = {:.2f}'.format(q, rho_at_cutoff))
    ls = '--' if q == 0.5 else '-.' 
    ax2.axhline(rho_at_cutoff, c='k', ls=ls)
    
    
    
ax0.set_title(r'simulated $R^2_{c,b}$')
ax1.set_title(r'empirical $R^2_{c,b}$')
ax2.set_title(r'inferred $\rho_{signal}$')

plt.savefig('PN_PC2_rho_signal.pdf', bbox_inches='tight')
plt.show()




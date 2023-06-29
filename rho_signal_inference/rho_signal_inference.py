# -*- coding: utf-8 -*-
"""
Rho inference
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd

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
############### Constants
##############################################################################

N_exps = 10000    
my_rhos = np.linspace(0, 1, 101)

RSQUARED_CALCIUM = 0.77
RSQUARED_BRP = 0.75
RSQUARED_OCT_AIR_3h = 0.27
RSQUARED_OCT_MCH_3h = 0.12



##############################################################################
############### Model information
##############################################################################


### OCT vs AIR

# ORNs
exp_name = 'OCT vs AIR' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'ORN PC 1'
model_RSQUARED = 0.23
n_flies = 30

exp_name = 'OCT vs AIR' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'ORN interpreted PC 1'
model_RSQUARED = 0.25
n_flies = 30

# PNs
exp_name = 'OCT vs AIR' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'PN PC 1'
model_RSQUARED = 0.084
n_flies = 53

exp_name = 'OCT vs AIR' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'PN interpreted PC 1'
model_RSQUARED = 0.098
n_flies = 53


### OCT vs MCH

# ORNs
exp_name = 'OCT vs MCH' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'ORN PC 1'
model_RSQUARED = 0.031
n_flies = 35

# ORN pre-synaptic density
exp_name = 'OCT vs MCH' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'brp' # 'brp' or 'calcium'
model_predictor_name = 'ORN Brp-Short PC 2'
model_RSQUARED = 0.088
n_flies = 53

# PNs
exp_name = 'OCT vs MCH' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'PN PC 2'
model_RSQUARED = 0.20
n_flies = 69

exp_name = 'OCT vs MCH' # 'OCT vs AIR' or 'OCT vs MCH'
predictor_type = 'calcium' # 'brp' or 'calcium'
model_predictor_name = 'PN interpreted PC 2'
model_RSQUARED = 0.12
n_flies = 69

##############################################################################
############### Run simulations
##############################################################################

if predictor_type == 'brp':
    r_predictor = RSQUARED_BRP**(1/4)
elif predictor_type == 'calcium':
    r_predictor = RSQUARED_CALCIUM**(1/4)

if exp_name == 'OCT vs AIR':
    r_b = RSQUARED_OCT_AIR_3h**(1/4)
elif exp_name == 'OCT vs MCH':
    r_b = RSQUARED_OCT_MCH_3h**(1/4)

predictor_color = d_color['ORN'] if 'ORN' in model_predictor_name else d_color['PN']

save_tag = (predictor_type + ' ' + model_predictor_name).replace(' ', '-')

predictor_rsquared_string = r'$R^2_{c,b}$' if predictor_type == 'calcium' else r'$R^2_{brp,b}$'

N_rhos = len(my_rhos)

calculated_rpredictorb2s = np.zeros((N_exps, N_rhos))

check_predictor_rsquareds = []
check_behav_rsquareds = []

for ii in range(N_exps):
    for rr in range(N_rhos):
        
        cur_rho = my_rhos[rr]
        
        r_squared_predictorb, x_predictor_latent, x_behav_latent, \
            x_predictor_obs1, x_predictor_obs2, x_behav_obs1, x_behav_obs2 \
                = run_forward_analysis(cur_rho, r_predictor, r_b, n_flies)
    
        check_predictor_rsquareds.append(np.corrcoef(x_predictor_obs1, x_predictor_obs2)[0,1]**2)
        check_behav_rsquareds.append(np.corrcoef(x_behav_obs1, x_behav_obs2)[0,1]**2)
        
        calculated_rpredictorb2s[ii, rr] = r_squared_predictorb
    

##############################################################################
############### Visualize outputs
##############################################################################

# check that R^2 between observed calcium/behavior matches empirical R^2
fig, axs = plt.subplots(1, 2, figsize=(10,5), sharey=True)
axs[0].hist(check_predictor_rsquareds, color=predictor_color, bins=50, alpha=0.5)
axs[0].axvline(r_predictor**4, c='k', ls='--')
axs[0].set_title('{} repeatability check'.format(predictor_type))
axs[0].set_xlabel(predictor_rsquared_string)
axs[0].set_ylabel('# simulations')
axs[1].hist(check_behav_rsquareds, color=predictor_color, bins=50, alpha=0.5, label='simulations')
axs[1].axvline(r_b**4, c='k', ls='--', label='empirical')
axs[1].set_title('behavior repeatability check')
axs[1].set_xlabel('behavior $R^2$')
axs[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0)
plt.show()


model_rsquared_slice = model_RSQUARED + 0.2 * model_RSQUARED * np.array([-1, 1])


# plot R^2_cb R^2_brpb vs. rho
fig, axs = plt.subplots(1, 1, figsize=(6, 5), sharex=True); axs = [axs]
axs[0].axvline(model_RSQUARED, color=predictor_color, ls='--')
cmap = colormaps['binary']
sim_res = [calculated_rpredictorb2s]
plot_maxx = 0
for rri in [0]:
    cur_ax = axs[rri]

    cur_calculated_rpredictorb2s = sim_res[rri]
    
    cur_median_rpredictorb2s = np.median(cur_calculated_rpredictorb2s, 0)
    for i in range(1, 10):
        cur_pctile = i/10 - 1/20
        c_arg = 1 - np.abs(cur_pctile - 0.5 + 0.05)*2
        my_clr = cmap(c_arg)
        
        cur_pctile_rsquared_maxes = np.quantile(cur_calculated_rpredictorb2s, cur_pctile, axis=0)
        next_pctile_rsquared_maxes = np.quantile(cur_calculated_rpredictorb2s, cur_pctile + 0.1, axis=0)        
        
        if i < 10:
            cur_ax.fill_betweenx(my_rhos, 
                             x1 = cur_pctile_rsquared_maxes, 
                             x2 = next_pctile_rsquared_maxes, 
                             color=my_clr, lw=0, alpha=0.55)
        
            maxval = np.max(np.concatenate((cur_pctile_rsquared_maxes, next_pctile_rsquared_maxes)))
            #print(i, maxval)
            plot_maxx = np.max([plot_maxx, maxval])
    
    cur_ax.plot(cur_median_rpredictorb2s, my_rhos, c='k')#, label='median', c='k')
    #cur_ax.legend(loc='upper left', bbox_to_anchor=(1.03, 1), borderaxespad=0)
for ax in axs:
    ax.set_ylim(0, 1)
    ax.set_xlim(0, plot_maxx)
axs[0].set_title('{}: {} ($R^2$ = {:.2f})\n'.format(exp_name, model_predictor_name, model_RSQUARED) + \
                 predictor_rsquared_string + ': {:.2f}'.format(r_predictor**4) + \
                 ', $R^2_{b,b}$: ' + \
                 r'{:.2f}, $N$: {:d}'.format(r_b**4, n_flies))
axs[0].set_ylabel(r'$\rho_{signal}$')
axs[0].set_xlabel(r'$R^2_{c,b}$' if predictor_type == 'calcium' else r'$R^2_{brp,b}$')

axs[0].fill_between(x=model_rsquared_slice, y1=1, y2=0, alpha=0.2, color=predictor_color)

plt.savefig('{}_{}_rho_vs_rsquared.png'.format(exp_name, save_tag), dpi=400, bbox_inches='tight')
plt.savefig('{}_{}_rho_vs_rsquared.pdf'.format(exp_name, save_tag))
plt.show()


#### marginal rho distribution
where_in_model_slice = (calculated_rpredictorb2s >= model_rsquared_slice[0]) & (calculated_rpredictorb2s <= model_rsquared_slice[1])
where_in_model_slice_marg = where_in_model_slice.sum(0)
where_in_model_slice_marg2 = where_in_model_slice_marg / np.trapz(where_in_model_slice_marg, x=my_rhos)

where_in_model_slice_p = where_in_model_slice_marg / np.sum(where_in_model_slice_marg)
where_in_model_slice_cump = np.cumsum(where_in_model_slice_p)

fig, axs = plt.subplots(1, 1, figsize=(6,5))#, sharex=True, gridspec_kw = {'height_ratios': [1, 0.4]})
axs.plot(my_rhos, where_in_model_slice_marg2, color=predictor_color, label=model_predictor_name)

qs = np.array([0.05, 0.5, 0.95])
rho_at_qs = np.zeros(len(qs))
for iq in range(len(qs)):
    lstyle = '-' if iq == 1 else '--'
    qi_rho_predictor = np.where(where_in_model_slice_cump <= qs[iq])[0][-1]
    rho_at_qs[iq] = my_rhos[qi_rho_predictor]
    axs.plot([my_rhos[qi_rho_predictor]]*2, 
                   [0, where_in_model_slice_marg2[qi_rho_predictor]], color=predictor_color, ls=lstyle)
    
plt.xlabel(r'$\rho_{signal}$')
axs.set_ylabel('marginal density')
axs.set_title('{}: {} ($R^2$ = {:.2f})'.format(exp_name, model_predictor_name, model_RSQUARED))

axs.set_ylim(0, axs.get_ylim()[1])
plt.xlim(0, 1)
plt.savefig('{}_{}_marginal_rho.png'.format(exp_name, save_tag), dpi=400, bbox_inches='tight')
plt.savefig('{}_{}_marginal_rho.pdf'.format(exp_name, save_tag))
plt.show()


##############################################################################
############### Visualize outputs
##############################################################################

row_index = ['experiment', 'predictor type', 'model predictor', 'save_tag', 
             'R^2 model', 'R^2 predictor', 'R^2 behav', 'num flies', 
             'rho_05', 'rho_50', 'rho_95']
row_info = [exp_name, predictor_type, model_predictor_name, save_tag, 
            model_RSQUARED, r_predictor**4, r_b**4, n_flies, 
            *rho_at_qs]
model_info = pd.Series(row_info, index=row_index)


# if df_model_rho_summary doesn't exist yet:
# df_model_rho_summary = pd.DataFrame(index=row_index)
# df_model_rho_summary.to_csv('df_model_rho_summary.csv')
df_model_rho_summary = pd.read_csv('df_model_rho_summary.csv')
df_model_rho_summary = df_model_rho_summary.append(model_info, ignore_index=True)
df_model_rho_summary.to_csv('df_model_rho_summary.csv', index=False)

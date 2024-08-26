'''
Code for plotting trial/lobe calcium responses in ORNs or PNs for all flies
(heatmaps and scatter plots)

Notes: calcium data loaded here was obtained by running 
ORNvsPN_analysis_ALLDATA/plotRandomSubsetPNandORN.m/evaluate_allData.m
and saving gh146rawdata and ocroRawData prior to NaN-infilling,
and also saving flyindices, flyindicesL, flyindicesR (PN fly/trial/lobe data),
and flyindicesO, flyindicesOL, flyindicesOR (ORN fly/trial/lobe data)
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr 

def set_font_sizes(SMALL_SIZE=14, MEDIUM_SIZE=16, LARGE_SIZE=20):
    '''
    Sets font size for matplotlib
    From: https://stackoverflow.com/a/39566040
    '''
    font = {'family':'sans-serif',
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

# set up ORN and PN colors
d_color = {'ORN': '#ffc951',
           'PN':  '#8621f1'}

gnames = ['DC2', 'DL5', 'DM1', 'DM2', 'DM3']

onames = ['air',
          '3-octanol',
          '1-hexanol',
          'ethyl lactate',
          'citronella',
          '2-heptanone',
          '1-pentanone',
          'ethanol',
          'geranyl acetate',
          'hexyl acetate',
          '4-methylcyclohexanol',
          'pentyl acetate',
          '1-butanol']

cnames = np.concatenate([[f'{x}_{y}' for y in onames] for x in gnames])


responsesPNs = pd.read_csv('gh146rawdata.csv', header=None)
flyindicesP = pd.read_csv('ALLPNflyindices.csv', header=None)
flyindicesPL = pd.read_csv('ALLPNflyindicesL.csv', header=None)
flyindicesPR = pd.read_csv('ALLPNflyindicesR.csv', header=None)

responsesORNs = pd.read_csv('orcorawdata.csv', header=None)
flyindicesO = pd.read_csv('ALLORNflyindices.csv', header=None)
flyindicesOL = pd.read_csv('ALLORNflyindicesL.csv', header=None)
flyindicesOR = pd.read_csv('ALLORNflyindicesR.csv', header=None)


df_responses_PN = responsesPNs.copy()
df_responses_PN.columns = cnames
df_responses_PN['fly'] = 0
df_responses_PN['lobe'] = 'X'
df_responses_PN['trial'] = -1


df_responses_ORN = responsesORNs.copy()
df_responses_ORN.columns = cnames
df_responses_ORN['fly'] = 0
df_responses_ORN['lobe'] = 'X'
df_responses_ORN['trial'] = -1


l_celltypes = ['PN', 'ORN']
l_df_responses = [df_responses_PN, df_responses_ORN]
l_indices = [flyindicesP, flyindicesO]
l_indicesL = [flyindicesPL, flyindicesOL]
l_indicesR = [flyindicesPR, flyindicesOR]


df_all_responses = []
df_all_responses_long = []

for i in range(2):

   
    df_responses = l_df_responses[i]
    flyindices = l_indices[i]
    flyindicesL = l_indicesL[i]
    flyindicesR = l_indicesR[i]
    cell_type = l_celltypes[i]
    

    
    df_fly_trials = []
    
    n_flies = len(flyindices)
    

    
    for i in range(n_flies):
        curL = flyindicesL.iloc[i].values
        curR = flyindicesR.iloc[i].values
        curind = flyindices.iloc[i].values
        
        curLexist = curL[~np.isnan(curL)]
        curRexist = curR[~np.isnan(curR)]
        curindexist = curind[~np.isnan(curind)]
        
        for j in range(len(curLexist)):
            val = curLexist[j]
            df_responses.loc[val-1, 'lobe'] = 'L'
            df_responses.loc[val-1, 'trial'] = j+1 
        
        for j in range(len(curRexist)):
            val = curRexist[j]
            df_responses.loc[val-1, 'lobe'] = 'R'
            df_responses.loc[val-1, 'trial'] = j+1 
            
        for val in curindexist:
            df_responses.loc[val-1, 'fly'] = i + 1 
            
        #behav_scores[i] = prefs[int(curind[0])-1]
        
    df_responses['cell_type'] = cell_type 
    
    df_trials_per_fly = (df_responses
                         .sort_values(['fly', 'lobe', 'trial'], ascending=[1,1,0])
                         .drop_duplicates(['fly', 'lobe'])
                         .pivot_table(index='fly', columns='lobe', values='trial')
                         .fillna(0)
                         .astype(int)
                         .assign(tot = lambda x: x['L'] + x['R'])
                         )
    
    df_responses_long = df_responses.melt(id_vars = ['fly', 'lobe', 'trial', 'cell_type'], value_name='response')
    df_responses_long[['glom', 'odor']] = df_responses_long['variable'].str.split('_', expand=True)      
    df_responses_long = df_responses_long.drop(columns='variable')    

    df_all_responses.append(df_responses)
    df_all_responses_long.append(df_responses_long)
    

df_all_responses_long = pd.concat(df_all_responses_long).reset_index(drop=True)



gh146tri = (df_all_responses_long
            .pipe(lambda x: x[x['cell_type'] == 'PN'])
            .dropna(subset='response')
            #.groupby(['fly', 'lobe', 'trial', 'glom', 'odor'])
            #.mean('response')
            .pivot_table(index=['fly', 'lobe', 'trial'], columns=['glom', 'odor'], values='response')
            )





lobedx = (df_all_responses_long
 .dropna(subset='response')
 )


df_lobes = (lobedx.loc[lobedx.lobe == 'R', ['cell_type', 'fly', 'glom', 'odor', 'response']]
           .merge(lobedx.loc[lobedx.lobe == 'L', ['cell_type', 'fly', 'glom', 'odor', 'response']],
                  on=['cell_type', 'fly', 'glom', 'odor'],
                  suffixes=['_R', '_L']
                  )
           )


df_trials = (df_all_responses_long.loc[df_all_responses_long.trial == 1]
         .merge(df_all_responses_long.loc[df_all_responses_long.trial == 2], 
                on=['cell_type', 'fly', 'lobe', 'glom', 'odor'],
                suffixes=['_1', '_2'])
         )


df_tri = df_trials.dropna(subset=['response_1', 'response_2'])



###################### SCATTER PLOT

tsize = 14
alph = 0.1 
fig, axs = plt.subplots(5, 5, 
                        height_ratios = [1, 1, 0.2, 1, 1], 
                        figsize=(12*1.2,12*1.2), sharex=True, sharey=True)

tlim = (-0.5, 4.5)
tt = np.arange(0, 4.1)
tt = [0, 4]
tt = [-0.5, 4.5]
tt = [-1, 5]
tlim = (-1, 5)

for iG in range(5):
    cur_glom = gnames[iG]
    
    ax = axs[0, iG]
    subdf = df_lobes[(df_lobes.cell_type == 'ORN') & (df_lobes.glom == cur_glom)]
    
    xx = subdf['response_R']
    yy = subdf['response_L']
    spear = spearmanr(xx, yy)
    
    sigcode = ''#'*' if np.log10(spear.pvalue) < -5 else ''
    ax.grid()
    ax.scatter(xx, yy,
               alpha=alph, lw=0, c=d_color['ORN'], zorder=10)
    ax.set_title(f'{cur_glom}\n' + \
                 r'$\rho = {:.2f}{}, n = {:.0f}$'.format(spear.statistic, sigcode, len(xx)), size=tsize)#r'$\rho = {:.2f}$ (n={})'.format(spear.statistic, len(xx)))
    #print(spear.pvalue)
    ax.set_xlim(tlim); ax.set_ylim(tlim)
    ax.set_xticks(tt); ax.set_yticks(tt)
    ax.xaxis.set_ticks_position('none'); ax.yaxis.set_ticks_position('none') 
    
        
    ax = axs[1, iG]
    subdf = df_lobes[(df_lobes.cell_type == 'PN') & (df_lobes.glom == cur_glom)]
    xx = subdf['response_R']
    yy = subdf['response_L']
    spear = spearmanr(xx, yy)
    
    sigcode = ''#'*' if np.log10(spear.pvalue) < -5 else ''
    ax.grid()
    ax.scatter(xx, yy,
               alpha=alph, lw=0, c=d_color['PN'], zorder=10)
    ax.set_title(#f'{cur_glom} ORNs\n' + 
                 r'$\rho = {:.2f}{}, n = {:.0f}$'.format(spear.statistic, sigcode, len(xx)), size=tsize)
    #print(spear.pvalue)
    ax.set_xlim(tlim); ax.set_ylim(tlim)
    ax.set_xticks(tt); ax.set_yticks(tt)
    ax.xaxis.set_ticks_position('none'); ax.yaxis.set_ticks_position('none') 
    
axs[1, 0].set_xlabel('right $\Delta$f/f')
axs[1, 0].set_ylabel('left $\Delta$f/f')

for ax in [axs[2, 0], axs[2, 1], axs[2, 2], axs[2, 3], axs[2, 4]]:
    ax.axis('off')


for iG in range(5):
    cur_glom = gnames[iG]
    
    ax = axs[3, iG]
    subdf = df_tri[(df_tri.cell_type == 'ORN') & (df_tri.glom == cur_glom)]
    xx = subdf['response_1']
    yy = subdf['response_2']
    
    spear = spearmanr(xx,yy)
    
    sigcode = ''#'*' if np.log10(spear.pvalue) < -5 else ''
    #ax.grid()
    ax.scatter(xx, yy, 
               alpha=alph, lw=0, c=d_color['ORN'], zorder=10)
    ax.set_title(f'{cur_glom}\n' + \
                 r'$\rho = {:.2f}{}, n = {:.0f}$'.format(spear.statistic, sigcode, len(xx)), size=tsize)
    ax.set_xlim(tlim); ax.set_ylim(tlim)
    ax.set_xticks(tt); ax.set_yticks(tt)
    ax.xaxis.set_ticks_position('none'); ax.yaxis.set_ticks_position('none') 
    
        
    ax = axs[4, iG]
    subdf = df_tri[(df_tri.cell_type == 'PN') & (df_tri.glom == cur_glom)]
    xx = subdf['response_1']
    yy = subdf['response_2']
    
    spear = spearmanr(subdf['response_1'], subdf['response_2'])
    
    sigcode = ''#'*' if np.log10(spear.pvalue) < -5 else ''
    #ax.grid()
    ax.scatter(xx, yy,
               alpha=alph, lw=0, c=d_color['PN'], zorder=10)
    ax.set_title(#f'{cur_glom} ORNs\n' + 
                 r'$\rho = {:.2f}{}, n = {:.0f}$'.format(spear.statistic, sigcode, len(xx)), size=tsize)
    #print(spear.pvalue)
    ax.set_xlim(tlim); ax.set_ylim(tlim)
    ax.set_xticks(tt); ax.set_yticks(tt)
    ax.xaxis.set_ticks_position('none'); ax.yaxis.set_ticks_position('none') 
    
axs[-1, 0].set_xlabel('trial 1 $\Delta$f/f')
axs[-1, 0].set_ylabel('trial 2 $\Delta$f/f')

plt.subplots_adjust(hspace=0.45)
#plt.tight_layout()

plt.savefig('scatter_lobes_trials.pdf')
plt.show()
    


########### HEATMAPS
 
# ORN: 65 flies, 208 trials
# PN: 122 flies, 406 trials
      



for ctype in ['ORN', 'PN']:
    df_responses_long = df_all_responses_long[df_all_responses_long.cell_type == ctype].reset_index(drop=True)

    ######## RECREATE CURRENT SUP FIG
    
    glom_odor_order = (df_responses_long 
                    .dropna(subset='response')
                    .drop_duplicates(subset=['glom', 'odor'])
                    .assign(glom_odor = lambda x: x['glom'] + '-' + x['odor'],
                            ordr = lambda x: np.arange(len(x)))
                    .loc[:, 'glom_odor']
                    .values
                    )
                    
    
    heatdf = (df_responses_long
              .dropna(subset='response')
              .assign(glom_odor = lambda x: x['glom'] + '-' + x['odor'],
                      ordr = lambda x: np.arange(len(x)))
              .groupby(['fly', 'glom_odor', 'ordr'])
              .mean('response')
              .reset_index()
              .pivot_table(index='glom_odor', columns='fly',  values='response')
              .loc[glom_odor_order, :]
              )
    
      
    
    fly_trial_order = (df_responses_long 
                    .dropna(subset='response')
                    .drop_duplicates(subset=['fly', 'lobe' , 'trial'])
                    .sort_values(['fly', 'lobe', 'trial'], ascending=[1,1,1])
                    .assign(fly_tri = lambda x: x['fly'].astype(str) + '_' + x['lobe'].astype(str) + '_' + x['trial'].astype(str))
                    .loc[:, 'fly_tri']
                    .values
                    )
    

    heatdf = (df_responses_long
              .dropna(subset='response')
              .assign(glom_odor = lambda x: x['glom'] + '-' + x['odor'],
                      fly_tri = lambda x: x['fly'].astype(str) + '_' + x['lobe'].astype(str) + '_' + x['trial'].astype(str),
                      ordr = lambda x: np.arange(len(x)))
              .pivot_table(index='glom_odor', columns='fly_tri', values='response')
              .loc[glom_odor_order, fly_trial_order]
              
              )
    
    
    
    plt.figure(figsize=(10,10))
    axh = sns.heatmap(heatdf, 
                      linewidth=0.5, linecolor='k',
                      cbar=False,
                      vmin=-0.53, vmax=4.3, cmap='hot')
    axh.collections[0].cmap.set_bad('0.5')
    plt.xlabel('individuals (each lobe/trial in separate columns)')
    plt.ylabel('glomerulus - odor')
    plt.xticks([])
    plt.yticks([])
    plt.show()
    
    
    plotdf = heatdf#.iloc[-10:, :30]
    ncol = plotdf.shape[1]
    
    fig, axs = plt.subplots(2, 1, figsize=(16/406*ncol,18), height_ratios=[75,1], sharex=True)
    
    axh = sns.heatmap(plotdf, 
                      linewidth=0.5,#.5, 
                      linecolor='k',
                      cbar=False,
                      vmin=-0.53, vmax=4.3, cmap='hot',
                      rasterized=False,
                      ax = axs[0])
    axh.collections[0].cmap.set_bad('0.5')
    
    
    
    margdf = pd.DataFrame(np.zeros((2, ncol)), columns=plotdf.columns)
    
    
    for i in range(ncol):
        cur_col = fly_trial_order[i]
        cur_fly, cur_lobe, cur_tri = re.findall('(\d+)_(\w)_(\d)', cur_col)[0]
        cur_fly = int(cur_fly)
        j = cur_fly % 2 
        c = 'C0' if cur_lobe == 'R' else 'C2'
        
        rect = Rectangle((i, j), 1, 1, 
                         linewidth=0, edgecolor='k', facecolor=c)
        axs[1].add_patch(rect)
        #margdf.iloc[j, i] = c
        
    lP, rP = axs[0].get_xlim()
    axs[0].set_xlim(lP - 0.5, rP + 0.5)
    
    axs[1].set_ylim(-0.1, 2.1)
    
    #axm = sns.heatmap(margdf, ax = axs[1], cbar=False, mask=margdf==0, 
    #                  cmap='Grays', vmin=0, vmax=1, annot=False)
    #axm.collections[0].cmap.set_bad('1')
    plt.xlabel('individuals (each lobe/trial in separate columns)')
    plt.ylabel('glomerulus - odor')
    #plt.xticks([])
    #plt.yticks([])
    
    axs[1].axis('off')
    
    
    for ax in axs:
        ax.xaxis.set_ticks_position('none'); ax.yaxis.set_ticks_position('none') 
        ax.set_yticks([]); ax.set_xticks([])
        ax.set_xlabel(''); ax.set_ylabel('')
        
    plt.subplots_adjust(hspace=0.01)
    plt.savefig(f'heatmap_flies_lobes_trials_{ctype}.pdf')
    plt.show()
    
    print(heatdf)
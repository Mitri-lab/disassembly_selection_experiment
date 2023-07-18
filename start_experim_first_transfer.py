import os
import sys
import pandas as pd
import selection_funcs_v010 as selfun

##### Define experimental parameters
spec_exclude = {'2':'122','122':('2','101'),'3':('1w','108'),'1w':'3',
                '101':'122','103':'114','108':'3','114':'103'}
specorder = [1,2,122,3,119,'1w',101,103,108,109,114]
Nspec  = 11  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 1   # Number of blanks (negative controls)
Ntop   = 10  # Number of communities to propagate to next round
Nspc   = 4

##### Collect first set of tubes
specpool = [str(s) for s in specorder]
params1 = {'Nspc':Nspc,'tubevec':range(1,Ntubes-Nblank+1),
           'specpool':specpool,'spec_exclude':spec_exclude}
params2 = {'Nspc':Nspc,'tubevec':range(Ntubes+1,2*Ntubes-Nblank+1),
           'specpool':specpool,'spec_exclude':spec_exclude}

df_treat1 = selfun.get_veryfirst_tubes(params1)
df_treat2 = selfun.get_veryfirst_tubes(params2)

#### First format: list composition by tube
outfile = sys.argv[1]
print('Writing first format, species by tubes, to '+outfile)
with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    df_treat1.astype('str').to_excel(writer,sheet_name='Treatment 1')
    df_treat2.astype('str').to_excel(writer,sheet_name='Treatment 2')

#### Second format: list tubes by species
outfile = sys.argv[2]
print('Writing second format, tubes by species, to '+outfile)
selfun.set_tubesbyspec((df_treat1.astype('str'),df_treat2.astype('str')),outfile)

#### Third output: draft of COD sheet
outfile = sys.argv[3]
print('Creating COD file template as '+outfile)
dfCODtmp = pd.DataFrame({'Tube':range(1,2*Ntubes+1),
                         'COD0':['']*(2*Ntubes),
                         'COD4':['']*(2*Ntubes)}).set_index('Tube')
with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    dfCODtmp.astype('str').to_excel(writer,sheet_name='COD input')

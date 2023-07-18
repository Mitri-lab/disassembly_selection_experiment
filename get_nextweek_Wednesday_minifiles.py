import sys
import string
import numpy as np
import pandas as pd
import selection_funcs_v010 as selfun
import openpyxl

from openpyxl.workbook import Workbook
from openpyxl.styles import PatternFill # Font, Fill,

##### Define experimental parameters
Nspec  = 11  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 1   # Number of blanks (negative controls)
Ntop   = 10  # Number of communities to propagate to next round
scorepow = 1.0
penpow = 1.0
N_spc = 4
pmigr = 0.75

##### Define paths for input files
#basedir  = '/Users/bvessman/gitfolder/Selection/data/'
if len(sys.argv)<2:
    print('Input error: You need to specify which input and output files to use. Example:')
    print('$ python3 get_nextweek_Wednesday ./path_to_input/presence_minifile.xlsx ./path_to_output/Species_by_tubes.xlsx ./path_to_output/Tubes_by_species.xlsx')
    quit()

#### Load dataset
coltup_Wed = ('Tube','Sp_Nr','Sp_Presence','COD_Day0','COD_Day4','Disassembly_Agar','Survival_Day4','Extin-Cont_Day4')

presfile = sys.argv[1]
# presfile = './Pilot2/20190815_4SC-ASE_Pilot2_Monday_dissassembly.xlsx'
print('Reading species presence from: ')
print(presfile)
df_pres_T1 = pd.read_excel(presfile, header=0, sheet_name='Treatment 1')
df_pres_T2 = pd.read_excel(presfile, header=0, sheet_name='Treatment 2')
df_surv    = pd.read_excel(presfile, header=0, sheet_name='Survival') # Sara added index_col=0 to avoid 'Unnamed:0' being included as a column

### Small cheat: pd.DataFrame.rename(columns={"old":"new"})
##  does nothing if the column 'old' is not present in the dataframe
##  For safety, if 'Unnamed: 0' then XYZ could be used
# if "Unnamed: 0" in df_pres_T1.columns:
df_pres_T1 = df_pres_T1.rename(columns={"Unnamed: 0": "Tube"})
df_pres_T2 = df_pres_T2.rename(columns={"Unnamed: 0": "Tube"})
df_surv = df_surv.rename(columns={"Unnamed: 0": "Tube"})

###### Step 1, Wednesday: Determine which tubes to disassemble ######
#### Collect presence matrix for disassembly
CODfile = sys.argv[2]
print('Reading COD from: ')
print(CODfile)
df_COD = pd.read_excel(CODfile, header=0)

outfile = presfile
print('Marking which species to collect from which plated community in: ')
print(outfile)
df_tuple = ((df_pres_T1.set_index('Tube'),df_pres_T2.set_index('Tube')),
             df_surv.set_index('Tube'),df_COD.set_index('Tube'))
selfun.colorExcel_specfromtubes_mf(df_tuple,[1,1-Ntubes],presfile,outfile)

###### Step 2, Wednesday: Determine which tubes to disassemble ######
df_master= selfun.get_master_df(df_tuple)
#print(df_master)

#### Define Survival, extinction etc to column 'Extin-Cont_Day4'
conditions = ('Sp_Presence==1 and Survival_Day4==1',
              'Sp_Presence==1 and Survival_Day4==0',
              'Sp_Presence==0 and Survival_Day4==1',
              'Sp_Presence==0 and Survival_Day4==0')
labeltup   = ('Survival','Exctinction','Contamination','Absent')
for cond, label in zip(conditions,labeltup):
    idces = df_master.query(cond).index
    for idx in idces:
        df_master.loc[idx,'Extin-Cont_Day4'] = label

#### In selection treatment: replicate communities based on COD
pardict_sel  = {'scorepow':scorepow,'penpow':penpow,'N_spc':N_spc,
                'tubevec':range(1,Ntubes-Nblank+1)}
pardict_ctrl = {'scorepow':scorepow,'penpow':penpow,'N_spc':N_spc,
                'tubevec':range(Ntubes+1,2*Ntubes-Nblank+1)}
df_sel  = df_master.loc[df_master.Tube<=Ntubes,:].copy()
df_ctrl = df_master.loc[df_master.Tube>Ntubes,:].copy()

dftubes_sel, dfpm_sel_nochange = selfun.selection_treatment_premigr_mf(df_sel,pardict_sel)
dfpm_ctrl_nochange             = selfun.control_treatment_premigr_mf(df_ctrl,pardict_ctrl)

#### Migrate species into a small proportion of these communities
specpool = tuple(np.unique([str(s) for s in df_master.Sp_Nr]))
spec_exclude = {'2':'122','122':('2','101'),'3':('1w','108'),'1w':'3',
                '101':'122','103':'114','108':'3','114':'103'}
pardict  = {'specpool':specpool,'N_spc':N_spc,'premix':pmigr,'specexcl':spec_exclude}

#### Finally - copy before-migration communities and migrate!
notpres_sel  = [1,]       # Dummy variable to allow the first loop
while len(notpres_sel)>0: # In the (unlikely!) case of species not present after migration
    df_premigr_sel = dfpm_sel_nochange.copy()
    notpres_sel    = selfun.migrate_members_new(df_premigr_sel,pardict)

notpres_ctrl = [1,]
while len(notpres_ctrl)>0:
    df_premigr_ctrl = dfpm_ctrl_nochange.copy()
    notpres_ctrl    = selfun.migrate_members_new(df_premigr_ctrl,pardict)

####### Then, define two output formats and write out results to xlsx
coltup = tuple(['Sp_'+str(s) for s in np.arange(1,N_spc+1)])

#### First format: list composition by tube
outfile = sys.argv[3]
print('Writing first format, species by tubes, to '+outfile)
with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    df_premigr_sel.loc[:,coltup].astype('str').to_excel(writer,sheet_name='Treatment 1')
    df_premigr_ctrl.loc[:,coltup].astype('str').to_excel(writer,sheet_name='Treatment 2')

#### Second format: list tubes by species
outfile = sys.argv[4]
print('Writing second format, tubes by species, to '+outfile)
selfun.set_tubesbyspec((df_premigr_sel,df_premigr_ctrl),outfile)

#### Third output: draft of COD sheet
if len(sys.argv)>5:
    outfile = sys.argv[5]
    print('Creating COD file template as '+outfile)
    dfCODtmp = pd.DataFrame({'Tube':range(1,2*Ntubes+1),
                             'COD0':['']*(2*Ntubes),
                             'COD4':['']*(2*Ntubes)}).set_index('Tube')
    with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
        dfCODtmp.astype('str').to_excel(writer,sheet_name='COD input')

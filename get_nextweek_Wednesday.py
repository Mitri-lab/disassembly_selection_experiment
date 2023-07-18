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
scorepow = 5.5
penpow = 1.0
N_spc = 4
pmigr = 0.5

##### Define paths for input files
#basedir  = '/Users/bvessman/gitfolder/Selection/data/'
if len(sys.argv)<4:
    print('Input error: You need to specify which input and output files to use. Example:')
    print('$ python3 get_nextweek_Wednesday ./path_to_input/input.xlsx ./path_to_output/Dissasembly_Monday.xlsx ./path_to_output/Species_by_tubes.xlsx ./path_to_output/Tubes_by_species.xlsx')
    quit()

#### Load dataset
coltup_Wed = ('Tube','Sp_Nr','Sp_Presence','COD_Day0','COD_Day4','Disassembly_Agar','Survival_Day4','Extin-Cont_Day4')
#datafile = '20190718_4SC-ASE_Pilot1.xlsx'

datafile = sys.argv[1]
print('Reading from '+datafile)
datadict_raw = pd.read_excel(datafile, header=0, sheet_name='Raw')

#### Split into selection and control
#### Soon, these will be treatment 1 and 2 to blind
df_sel  = datadict_raw.loc[datadict_raw.Selection_Treatment=='Group',coltup_Wed].copy()
df_ctrl = datadict_raw.loc[datadict_raw.Selection_Treatment=='Random',coltup_Wed].copy()
df_tuple = (df_sel,df_ctrl)

###### Step 1, Wednesday: Determine which tubes to disassemble ######
#### Collect presence matrix for disassembly
infile  = sys.argv[2] # './Table_disassembly_Monday.xlsx'
outfile = infile
print('Marking which species to collect from which plated community in '+infile)
selfun.colorExcel_specfromtubes(df_tuple,[1,1-Ntubes],infile,outfile)

###### Step 2, Wednesday: Determine which tubes to disassemble ######
#### Define Survival, extinction etc to column 'Extin-Cont_Day4'
conditions = ('Sp_Presence==1 and Survival_Day4==1',
              'Sp_Presence==1 and Survival_Day4==0',
              'Sp_Presence==0 and Survival_Day4==1',
              'Sp_Presence==0 and Survival_Day4==0')
labeltup   = ('Survival','Exctinction','Contamination','Absent')
for df in (df_sel,df_ctrl):
    for cond, label in zip(conditions,labeltup):
        idces = df.query(cond).index
        for idx in idces:
            df.loc[idx,'Extin-Cont_Day4'] = label

#### In selection treatment: replicate communities based on COD
pardict= {'scorepow':scorepow,'penpow':penpow,'N_spc':N_spc}
dftubes_sel, dfpm_sel_nochange = selfun.selection_treatment_premigr(df_sel,pardict)
dfpm_ctrl_nochange             = selfun.control_treatment_premigr(df_ctrl,pardict)

#### Migrate species into a small proportion of these communities
specpool = tuple(np.unique([str(s) for s in df_sel.Sp_Nr]))
spec_exclude = {'2':'122','122':('2','101'),'3':('1w','108'),'1w':'3',
                '101':'122','103':'114','108':'3','114':'103'}
pardict  = {'specpool':specpool,'N_spc':N_spc,'premix':pmigr,'specexcl':spec_exclude}

#### Finally - copy before-migration communities and migrate!
notpres_sel  = [1,]       # Dummy variable to allow the first loop
while len(notpres_sel)>0: # In the (unlikely!) case of species not present after migration
    df_premigr_sel = dfpm_sel_nochange.copy()
    notpres_sel  = selfun.migrate_members_new(df_premigr_sel,pardict)

notpres_ctrl = [1,]
while len(notpres_ctrl)>0:
    df_premigr_ctrl = dfpm_ctrl_nochange.copy()
    notpres_ctrl = selfun.migrate_members_new(df_premigr_ctrl,pardict)

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

###### Step 3, Wednesday: Define next week's Excel tables ######

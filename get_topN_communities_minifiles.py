import sys
import string
import numpy as np
import pandas as pd
import selection_funcs_v010 as selfun
import openpyxl

from pathlib import Path
from os.path import isfile, join
from openpyxl.workbook import Workbook
from openpyxl.styles import PatternFill # Font, Fill,

##### Define experimental parameters
from disassembly_parameters import *

##### Define columns and header strings
coltup_Wed = ('Tube','Sp_Nr','Sp_Presence','COD_Day0','COD_Day4',
              'Disassembly_Agar','Survival_Day4','Extin-Cont_Day4')

conditions = ('Sp_Presence==1 and Survival_Day4==1',
              'Sp_Presence==1 and Survival_Day4==0',
              'Sp_Presence==0 and Survival_Day4==1',
              'Sp_Presence==0 and Survival_Day4==0')
labeltup   = ('Survival','Exctinction','Contamination','Absent')

pardict_sel  = {'scorepow':scorepow,'penpow':penpow,'N_spc':N_spc,
                'tubevec':range(1,Ntubes-Nblank+1)}
pardict_ctrl = {'scorepow':scorepow,'penpow':penpow,'N_spc':N_spc,
                'tubevec':range(Ntubes+1,2*Ntubes-Nblank+1)}

###### Define output folder, file names and output format
outcols = ('Week','COD_Day0','COD_Day4','CODscore','Survival',
           'Sp_1','Sp_2','Sp_3','Sp_4','Score','Prob')

dfout_sel  = pd.DataFrame(columns=outcols)
dfout_ctrl = pd.DataFrame(columns=outcols)


#### Load dataset
for idxd, date in enumerate(datetup):
    print('Processing data from '+date)

    datadir = date+'_4SC-ASE_week'+str(idxd)
    datafile_pres = date+'_4SC-ASE_week'+str(idxd)+'_monday_disassembly.xlsx'
    datafile_COD  = date+'_4SC-ASE_week'+str(idxd)+'_COD.xlsx'

    presfile = join(basedir,datadir,datafile_pres)
    CODfile  = join(basedir,datadir,datafile_COD)

    if not Path(presfile).is_file():
        print('Error! Different name format for presence-absence file')
        print('Skipping date '+date)
        continue

    if not Path(CODfile).is_file():
        print('Error! Different name format for COD file')
        print('Skipping date '+date)
        continue

    #### Load presence-absence data
    # print('Reading species presence from: ')
    # print(presfile)
    df_pres_T1 = pd.read_excel(presfile, header=0, sheet_name='Treatment 1')#.query('Disassembly_Agar>0')
    df_pres_T2 = pd.read_excel(presfile, header=0, sheet_name='Treatment 2')#.query('Disassembly_Agar>0')
    df_surv    = pd.read_excel(presfile, header=0, sheet_name='Survival')

    ### Small cheat: pd.DataFrame.rename(columns={"old":"new"})
    ##  does nothing if the column 'old' is not present in the dataframe
    ##  For safety, if 'Unnamed: 0' then XYZ could be used
    # if "Unnamed: 0" in df_pres_T1.columns:
    df_pres_T1 = df_pres_T1.rename(columns={"Unnamed: 0": "Tube"})
    df_pres_T2 = df_pres_T2.rename(columns={"Unnamed: 0": "Tube"})
    df_surv = df_surv.rename(columns={"Unnamed: 0": "Tube"})

    #### Load also COD
    # print('Reading COD from: ')
    # print(CODfile)
    df_COD = pd.read_excel(CODfile, header=0)

    # print(df_pres_T1)
    ###### Get Master data file with score and extinction
    df_tuple = ((df_pres_T1.set_index('Tube'),df_pres_T2.set_index('Tube')),
                 df_surv.set_index('Tube'),df_COD.set_index('Tube'))
    df_master = selfun.get_master_df(df_tuple)

    #### Define Survival, extinction etc to column 'Extin-Cont_Day4'
    for cond, label in zip(conditions,labeltup):
        idces = df_master.query(cond).index
        for idx in idces:
            df_master.loc[idx,'Extin-Cont_Day4'] = label

    #### In selection treatment: replicate communities based on COD
    df_sel  = df_master.loc[df_master.Tube<=Ntubes,:].copy()
    df_ctrl = df_master.loc[df_master.Tube>Ntubes,:].copy()

    dftubes_sel,  dfpm_sel_nochange  = selfun.selection_treatment_premigr_mf(df_sel,pardict_sel)
    dftubes_ctrl, dfpm_ctrl_nochange = selfun.selection_treatment_premigr_mf(df_ctrl,pardict_ctrl)

    dftubes_sel.loc[:,'Week']  = np.repeat([idxd,],len(dftubes_sel.index))
    dftubes_ctrl.loc[:,'Week'] = np.repeat([idxd,],len(dftubes_ctrl.index))

    dfout_sel  = dfout_sel.append(dftubes_sel,sort=False)
    dfout_ctrl = dfout_ctrl.append(dftubes_ctrl,sort=False)

###### Set 'week' as index
dfout_sel.loc[:,'Tube'] = dfout_sel.index
dfout_sel = dfout_sel.loc[:,printcols].set_index('Week')

dfout_ctrl.loc[:,'Tube'] = dfout_ctrl.index
dfout_ctrl = dfout_ctrl.loc[:,printcols].set_index('Week')

###### Print to Excel file
# outfile = sys.argv[3]
outpath = join(basedir,outfile)
print('Printing to file ')
print(outpath)

with pd.ExcelWriter(outpath,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    dfout_sel.astype('str').to_excel(writer,sheet_name='Treatment 1')
    dfout_ctrl.astype('str').to_excel(writer,sheet_name='Treatment 2')

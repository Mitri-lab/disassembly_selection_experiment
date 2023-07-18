import os
import sys
import numpy as np
import pandas as pd
import selection_funcs_v010 as selfun

##### Define experimental parameters
# specorder = ['1','2','122','3','119','1w','101','103','108','109','114']
specorder = [1,2,122,3,119,'1w',101,103,108,109,114]
Nspec  = 11  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 2   # Number of blanks (negative controls)
Ntop   = 10  # Number of communities to propagate to next round

##### Define paths for input files
#basedir  = '/Users/bvessman/gitfolder/Selection/data/'
if len(sys.argv)<3:
    print('Input error: You need to specify which input and output files to use. ')
    print('$ python3 get_disassembly_Monday ./path_to_input/input.xlsx ./path_to_output/output.xlsx')
    quit()

#### Load dataset
# collabels = ('D:H,J:L')
collabels = ('A:R')
#datafile = '20190718_4SC-ASE_Pilot1.xlsx'
#datadict_raw = pd.read_excel(join(basedir,datafile), header=0, sheet_name='Raw', usecols=collabels)
datafile = sys.argv[1]
print('Reading from '+datafile)
datadict_raw = pd.read_excel(datafile, header=0, sheet_name='Raw', usecols=collabels)

#### Split into selection and control
coltup_Mon = ('Tube','Sp_Nr','Sp_Presence','COD_Day0','COD_Day4','Disassembly_Agar')

#### Soon, these will be treatment 1 and 2 to blind
df_sel  = datadict_raw.loc[datadict_raw.Selection_Treatment=='Group',coltup_Mon].copy()
df_ctrl = datadict_raw.loc[datadict_raw.Selection_Treatment=='Random',coltup_Mon].copy()

######### Step 1, Monday: Determine which tubes to disassemble #########
#### Collect presence matrix for disassembly
pardict = {'Nspec':Nspec,'Ntop':Ntop}
dfpiv_sel  = selfun.make_disassembly_dfs(df_sel,pardict)

#### For control treatment: disassemble randomly
df_ctrl_tmp = df_ctrl.copy()
df_ctrl_tmp.loc[:,'COD_Day4'] = np.random.permutation(df_ctrl_tmp.COD_Day4)
dfpiv_ctrl = selfun.make_disassembly_dfs(df_ctrl_tmp,pardict)

##### Write out a suggestion for Wednesday presence
dfsel_disass  = dfpiv_sel.loc[dfpiv_sel.Disassembly_Agar>0,specorder]
dfctrl_disass = dfpiv_ctrl.loc[dfpiv_ctrl.Disassembly_Agar>0,specorder]
df_survival = dfsel_disass.append(dfctrl_disass).replace(1,0).replace(0,'')

#### Write out results to xlsx
# outfile = './disassembly3_190804.xlsx'
outfile = sys.argv[2]
with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    dfpiv_sel.to_excel(writer,sheet_name='Treatment 1')
    dfpiv_ctrl.to_excel(writer,sheet_name='Treatment 2')
    df_survival.to_excel(writer,sheet_name='Survival')

print('Saving disassembly file as:')
print(outfile)
selfun.cellcolorExcel_disassembly(outfile,(dfpiv_sel,dfpiv_ctrl),('E2EFDA','E25EA6'),[1,1-Ntubes])

#### If additional arguments pointing to CFU and COD were passed, append info to data
if len(sys.argv)>3 & len(sys.argv)<5:
    #### Read input
    infile_CFU = sys.argv[3] # infile_CFU = './CFU_input.csv'
    infile_COD = sys.argv[4] # infile_COD = './COD_input.csv'
    df_CFU = pd.read_csv(infile_CFU,sep=',',header=0)
    df_COD = pd.read_csv(infile_COD,sep=',',header=0)

    print('Reading CFU from:')
    print(infile_CFU)
    print('Reading COD from:')
    print(infile_COD)
    print('Updating file {0} with CFU, COD.'.format(datafile))
    # print(infile_COD)
    rw = selfun.readwrite_CFU_COD(datadict_raw,(df_CFU,df_COD),(dfpiv_sel,dfpiv_ctrl),datafile)

import os
import sys
import numpy as np
import pandas as pd
import selection_funcs_v010 as selfun

##### Define experimental parameters
# specorder = ['1','2','122','3','119','1w','101','103','108','109','114']
specorder = [1,2,122,3,119,'1w',101,103,108,109,114]
Nspec  = len(specorder)  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 1
Ntop   = 10  # Number of communities to propagate to next round
N_spc  = 4

##### Define paths for input files
#basedir  = '/Users/bvessman/gitfolder/Selection/data/'
if len(sys.argv)<4:
    print('Input error: You need to specify which input and output files to use. ')
    print('$ python3 get_disassembly_Monday_minifiles ./path_to_input/COD_input.xlsx ./path_to_input/Assembly_last_Thursday.xlsx ./path_to_output/Disassembly_Monday.xlsx')
    print('See example files in same folder as script. ')
    quit()

#### Load dataset
CODfile = sys.argv[1]
print('Reading COD from: '+CODfile)
df_CODin      = pd.read_excel(CODfile, header=0, sheet_name='COD input')

#### Treatment 1 and 2 in tubes 1-30 and 31-60, respectively
commfile = sys.argv[2]
print('Reading community composition from: '+commfile)
dfpres_in_T1  = pd.read_excel(commfile, header=0, sheet_name='Treatment 1')
dfpres_in_T2  = pd.read_excel(commfile, header=0, sheet_name='Treatment 2')

#### If needed, rename unnamed tube column
dfpres_in_T1 = dfpres_in_T1.rename(columns={"Unnamed: 0": "Tube"})
dfpres_in_T2 = dfpres_in_T2.rename(columns={"Unnamed: 0": "Tube"})

#### Set tube column as index
dfpres_in_T1 = dfpres_in_T1.set_index('Tube')
dfpres_in_T2 = dfpres_in_T2.set_index('Tube')

#### Treatment 1 collected as presence and COD
df_T1 = df_CODin.set_index('Tube').loc[dfpres_in_T1.index,].copy()
dfpres_in_T1.loc[:,'COD_Day4'] = df_T1.COD4

#### Treatment 2 collected as presence and !! RANDOMLY ORDERED !! COD
df_T2 = df_CODin.set_index('Tube').loc[dfpres_in_T2.index,].copy()
COD_permut = np.random.permutation(df_T2.COD4)
dfpres_in_T2.loc[:,'COD_Day4'] = pd.Series(index=df_T2.index,data=COD_permut)

######### Step 1, Monday: Determine which tubes to disassemble #########
#### Collect presence matrix for disassembly
pardict = {'Nspec':Nspec,'Ntop':Ntop,'Ntubes':Ntubes,'N_spc':N_spc,
           'Nblanks':Nblank,'speclist':[str(s) for s in specorder]}
dfpiv_T1 = selfun.make_disassembly_dfs_minifile(dfpres_in_T1,pardict)
dfpiv_T2 = selfun.make_disassembly_dfs_minifile(dfpres_in_T2,pardict)

##### Write out a suggestion for Wednesday presence
speclist = [str(s) for s in specorder]
dfT1_disass = dfpiv_T1.loc[dfpiv_T1.Disassembly_Agar>0,speclist]
dfT2_disass = dfpiv_T2.loc[dfpiv_T2.Disassembly_Agar>0,speclist]
df_survival = dfT1_disass.append(dfT2_disass).replace(1,0).replace(0,'')

#### Write out results to xlsx
# outfile = './disassembly3_190804.xlsx'
outfile = sys.argv[3]
with pd.ExcelWriter(outfile,engine='openpyxl',date_format='YYYYMMDD', mode='w') as writer:
    dfpiv_T1.to_excel(writer,sheet_name='Treatment 1')
    dfpiv_T2.to_excel(writer,sheet_name='Treatment 2')
    df_survival.to_excel(writer,sheet_name='Survival')

print('Saving disassembly file as:')
print(outfile)
selfun.cellcolorExcel_disassembly(outfile,(dfpiv_T1,dfpiv_T2),('E2EFDA','E25EA6'),[1,1-Ntubes])

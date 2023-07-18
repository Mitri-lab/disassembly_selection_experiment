import sys
import pandas as pd
import selection_funcs_v010 as selfun

##### Define experimental parameters
Nspec  = 11  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 2   # Number of blanks (negative controls)
Ntop   = 10  # Number of communities to propagate to next round

##### Define paths for input files
#basedir  = '/Users/bvessman/gitfolder/Selection/data/'
if len(sys.argv)<3:
    print('Input error: You need to specify which input and output files to use. ')
    print('$ python3 get_disassembly_Monday ./path_to_input/input.xlsx ./path_to_disassembly/disassemble.xlsx ./path_to_output/output.xlsx')
    quit()

#### Load dataset
collabels = ('D:H,J:M')
coltup_Wed = ('Tube','Sp_Nr','Sp_Presence','COD_Day0','COD_Day4','Disassembly_Agar','Survival_Day4','Extin-Cont_Day4')
#datafile = '20190718_4SC-ASE_Pilot1.xlsx'
#datadict_raw = pd.read_excel(join(basedir,datafile), header=0, sheet_name='Raw', usecols=collabels)

datafile = sys.argv[1]
print('Reading from '+datafile)
datadict_raw = pd.read_excel(datafile, header=0, sheet_name='Raw', usecols=collabels)

# print(datadict_raw.head(15))
#### Split into selection and control
#### Soon, these will be treatment 1 and 2 to blind
df_sel  = datadict_raw.loc[datadict_raw.Selection_Treatment=='Group',coltup_Wed].copy()
df_ctrl = datadict_raw.loc[datadict_raw.Selection_Treatment=='Random',coltup_Wed].copy()
df_tuple = (df_sel,df_ctrl)

######### Step 1, Wednesday: Determine which tubes to disassemble #########
#### Collect presence matrix for disassembly

infile  = sys.argv[2] # './test_FAS_AW_Mon.xlsx'
outfile = sys.argv[3] # './propagatethis_190806.xlsx'
selfun.colorExcel_specfromtubes(df_tuple,[1,1-Ntubes],infile,outfile)

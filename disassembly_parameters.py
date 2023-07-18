##### Define
### disassembly_parameters.py

Nspec  = 11  # Number of species (before evolution)
Ntubes = 30  # Number of tubes per experiment condition
Nblank = 1   # Number of blanks (negative controls)
Ntop   = 10  # Number of communities to propagate to next round
scorepow = 1.0
penpow   = 1.0
N_spc = 4    # Number of species per community
pmigr = 0.7

#### Some species cannot be distinguished on regular plates,
##   Hence, the following combinations are avoided
spec_exclude = {'2':'122',
                '122':('2','101'),
                '3':('1w','108'),
                '1w':'3',
                '101':'122',
                '103':'114',
                '108':'3',
                '114':'103'}

###################### FOR
### Define paths for input files,
basedir = '/Users/bvessman/gitfolder/Selection/4SC-ASE_AllWeeksArchive/'

datetup = ('20191010','20191017','20191024','20191114','20191121',
           '20191128','20191212','20200109','20200116','20200123',
           '20200130','20200206','20200213','20200220','20200227',
           '20200305','20200507','20200514','20200521')

###### Define output files
printcols = ('Week', 'Tube', 'Score', 'Sp_1', 'Sp_2', 'Sp_3', 'Sp_4', 
             'CODscore', 'Survival')
outfile = '4SC-ASE_community_scores_Week0to19_20200728.xlsx'

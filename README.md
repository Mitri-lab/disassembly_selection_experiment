# disassembly_selection_experiment
Scripts used in a community breeding experiment to disassemble bacterial communities over several weeks.

---------------- README for disassembly scripts ---------------
------ These include 'get_disassembly_Monday_minifiles.py'
------               'get_nextweek_Wednesday_minifiles.py'
------ The supporting functions are found in 'selection_funcs_v010.py'
---- by Björn Vessman (Bjorn.Vessman@unil.ch)
v.01.0 on 2019-08-05, basic form of each script
          2019-08-12, added functions to create a new data file for 'next' week and
                      modify it with new CFU, COD
          2019-08-26, Minifiles version for the Monday, Wednesday scripts.
          2019-10-04, Method updates



------ INDEX
--- 0: What are these files?
--- 1: Software pre-requisites
--- 2: Format of the input data
--- 3: Schedule - when and how to run the scripts
---- 3x: Step 00 - start_experim_first_transfer.py
---- 3a: Step 1  - get_disassembly_Monday.py
---- 3a1: Step 1 - get_disassembly_Monday_minifiles.py
---- 3b: Step 2  - get_nextweek_Wednesday.py
---- 3b1: Step2  - get_nextweek_Wednesday_minifiles.py
--- 98: TODO
--- 99: References


--- 0: What are these files?
The selection experiment by Flor Arias-Sanchez, Sara Mitri, Alice Wallef Géraldine Alberti, Samuele Testa and Bjorn Vessman is performed by means of a method to artificially select and 'breed' bacterial communities for degradation of Metal-Working Fluids (experimental system presented in [1]). The selection method takes into account the COD scores measured in the previous week and suggests a community composition for next week, in order to increase substrate degradation and community cooperation. The scores are assumed to be found in an Excel file of a particular structure (see "3: Input format" below) and the method is implemented in python for easy automation and randomisation of the communities.

In this README, there is the occasional suggestion to run scripts. In these cases, the terminal is denoted as "$", so "$ python helloworld.py" means to open a terminal window and write "python helloworld.py" in order to run the script "helloworld.py" with the python interpreter.

The folder contains a few python files, namely
- selection_funcs_v010.py, a collection of functions used by the scripts.
  This file is not a script in itself and cannot be run as it is written now.
- get_disassembly_Monday.py (deprecated)
  The script reads an Excel file with data on community composition and COD, and
  writes an Excel file with a suggested disassembly - how to recapture each species
  from the collection of tubes. The script is run after the COD measurements are noted,
  typically on a Monday. If wanted, the script can alter the main data file and update it
  with the inoculum CFU (day 0) and COD for day 0 and day 4.
- get_nextweek_Wednesday.py (deprecated)
  The script reads the main Excel data file and the Excel file output by the Monday script. Based on these, we mark which tube to sample each species from based on their survival. Each species should be sampled from the highest-scoring community in which it survived. The output is a coloured version of the Monday output (either in a new, separate file or in the same file) which marks where to sample each species.

--- 1: Software pre-requisites
The script is written for python 3.7 and is supposed to be run in a UNIX environment (a bash terminal or similar). Linux and Mac OS support UNIX shell terminals by default, but if your computer uses MS Windows, I would suggest dual-booting with a linux-based OS. The IT support of your institution or your go-to computer person will be happy to help. =)

If you are not sure which version of python that you are using, run
>> python --version
and just to be safe, use the alias 'python3' as the interpreter.

The script requires the following python libraries
- sys, os (basic system functions to interact with the computer, always available)
- scipy and numpy (for SCIentific and NUMerical computations)
- pandas (for manipulation of data tables)
- openpyxl (for read/write of MS Excel files)
---- Possibly, openpyxl needs a supporting library xlrd (Pronounced Excel-read and written to do exactly that)

If any of these libraries (aka 'modules') are missing in your python installation, you would see an error message like
>> ModuleNotFoundError: No module named 'EXAMPLE_MODULE_NAME'
... and you can in that case install them by running
>> $ pip3 install openpyxl

If the problem is on an even higher level - that you do not have access to the package manager pip, first run
>> $ python -m ensurepip

For reference, there is online documentation of each module. Use your search engine of choice and search for 'python pandas documentation' or similar.


--- 2: Format of the input data
---- Big data file '20190718_4SC-ASE_Pilot1.xlsx' or similar
The input file should be an excel file with Ntubes x Nspecies rows (in our experiment, we expect 2x30 tubes and 11 species, giving 660 rows in total). The columns should include, in the following order,
-Assembly_Day          Day of assembly of community, a Thursday with a YYYYMMDD date.
-Experiment            Experiment label '4SC-ASE' (four-species communities, artificial selection experiment)
-Round_nr             'Pilot2', 'Week45' or similar
-Selection_Treatment  'Group', 'Control' or 'Treatment1', 'Treatment2' if randomised
-Tube                  1, 2, ..., Ntubes
-Sp_Nr                 Species labels (not necessarily numbers) 1, 2, ..., 119, 1w
-Sp_Presence           Is the corresponding species present (1) or not (0)?
-COD_Day0              COD measurement on day 0
-log10CFUml_Day0       CFU estimation, day 0
-COD_Day4              COD measurement on day 4
-Disassembly_Agar      Should this tube be plated for community disassembly? (1 else 0)
-Survival_Day4         Is the species present (1) at day 4 or not (0)?
-Extin-Cont_Day4       Combining 'Sp_Presence' and 'Survival_Day4', we can find extinction events and contamination.
-Dilution_Day4         Lowest dilution at which colonies were found
-Droplet_Day4          Size of pipetted droplet
-CFU_Day4              Colonies counted at the dilution 'Dilution_Day4'
-log10CFUml_Day4       Estimate CFU based on 'CFU_Day4', 'Droplet_Day4' and 'Dilution_Day4'
-Selected4NextRound    Marker (1 or 0) from which tube to propagate colonies to next round
---- Input_CFUCOD.csv

--- 3: Schedule - when and to run the scripts
The selection method workflow is implemented in the python scripts, which are supposed to run when the corresponding data is available. The experimental procedure starts by inoculating the relevant communities on a Thursday (the communities is assumed to be known from last week's results), measuring the initial COD and noting which species is present in which community. On a Monday, COD measurements are taken and noted in the data file which is then read by 'get_disassembly_Monday.py'. On the following Wednesday, the disassembly tells us the survival of species in each plated community, and where we should sample each species. Finally, the COD measurements and community survival is used to propose new community compositions for next week.

Having installed all the relevant modules, the scripts are run as follows.

---- 3x: Step 0 - start_experim_first_transfer.py
At some point you would want to start, or start over, the experiment and on that day you need a suggestion for which species goes into which tube. Run 'start_experim_first_transfer.py' as follows
>> $ python3 start_experim_first_transfer.py ../path_to_output/Species_by_tubes_w1.xlsx ../path_to_output/Tubes_by_species_w1.xlsx

---- 3a1: Step 1v1 - get_disassembly_Monday.py
>> $ python3 get_disassembly_Monday.py ./path_to_input_data/2019mmdd_4SC-ASE_WeekN.xlsx ./path_to_disassembly/2019mmdd_4SCASE_WeekN_Disassembly_Monday.xlsx

The 'Monday' script runs with two arguments as above, and you need to provide paths to both the input and output files. The output file does not need to exist yet, and can be created by the script in the provided filename. For example, we have input in '20190718_4SC-ASE_Pilot1.xlsx' and define the output in '20190722_4SCASE_Pilot1_Disassembly_Monday.xlsx'

If you wish to update the main data file with the initial CFU, and COD measurements for day 0 and 4, define two .csv files (regular text files with values separated by commas, like the provided CFU_input_example.csv and COD_input_example.csv) and run the script as
>> $ python3 get_disassembly_Monday.py ../Data/20190718_4SC-ASE_Pilot1.xlsx ../Outs/20190722_4SCASE_Pilot1_Disassembly_Monday.xlsx ../Data/CFU_input_example.csv ../Data/COD_input_example.csv

---- 3a1: Step 1 - get_disassembly_Monday_minifiles.py
>> $ python3 get_disassembly_Monday_minifiles.py ./INPUT/COD_day0_day4_example.xlsx ./INPUT/Species_by_tube_example.xlsx ../OUTPUT/Dissassembly_minifiles_example.xlsx

The 'Monday' script runs with three arguments as above, and you need to provide paths to both the COD input and species presence, as well as the output files. The output file does not need to exist yet, and can be created by the script in the provided filename. For example, we have COD input in 'COD_day0_day4_example.xlsx' (see also example Excel tables in the script folder) and presence (defined by the Wednesday script last week) in 'Species_by_tube_example.xlsx' and define the output in 'Dissassembly_minifiles_example.xlsx'. Except for the input format, the scripts do the follow the same procedure as the scripts that worked with the Excel 'master' file.

---- 3b: Step 2 - get_nextweek_Wednesday.py
The Wednesday script is run with _four_ arguments - the master data file, the disassembly file produced by the Monday script, and a path to the file where you want to store the which-species-from-which-tube table. For example, we have input in '20190815_4SC-ASE_Pilot2_MASTER.xlsx', the previously constructed disassembly file 'Monday_dissassembly.xlsx' and define the two output formats in './Output/Species_by_tubes.xlsx', './Output/Tubes_by_species.xlsx'
>> $ python3 get_nextweek_Wednesday.py ./Input/20190815_4SC-ASE_Pilot2_MASTER.xlsx ./Input/Monday_dissassembly.xlsx ./Output/Species_by_tubes.xlsx ./Output/Tubes_by_species.xlsx

The script reads COD and community composition from the 'master' file and marks which species should be plated from which tube. Given COD and survival data, we compute next week's community composition and print it to Excel files. 

---- 3b1: Step 2 - get_nextweek_Wednesday_minifiles.py
>> $ python3 get_nextweek_Wednesday_minifiles.py ./Input/Monday_dissassembly.xlsx ./Minifiles/COD_day0_day4.xlsx ./Output/Species_by_tube.xlsx ./Output/Tubes_by_spec.xlsx
>> $ python3 get_nextweek_Wednesday_minifiles.py ./INPUT/Monday_disassembly.xlsx ./INPUT/COD_week1.xlsx ./OUTPUT/Species_by_Tubes.xlsx ./OUTPUT/Tubes_by_species.xlsx ./Next_week/INPUT/COD_week2.xlsx


A minifile version of 'get_nextweek_Wednesday.py', using the COD minifile (also used on Monday, see section 3a1) and disassembly file with three sheets: Treatment 1, Treatment 2 and Survival. Each sheet contains 1s and 0s marking whether the species is present or absent in that tube.



--- 98: TODO
There are a few hacks and temporary solutions since the first pilot experiment data does not exactly follow the assumptions.
-- selection_funcs_v010.py
- set_tubesbyspec(dfs,outfile)
cell widths in output Excel file



--- 99: References
[1] Piccardi, Vessman, Mitri (2019) 'Toxicity drives facilitation between four bacterial species' PNAS

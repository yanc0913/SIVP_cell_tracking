# SIVP_cell_tracking
These two scripts are used to track and analyse endothelial cell (EC) migration during SIVP development.

Cell migratory statistics were obtained from TrackMate (Tinevez et al., 2017), Fiji, including Track statistics, Links in tracks statistics (referred as **links file**), and Spots in tracks statistics (referred as **spots file**). 

- Every individual EC traced was represented as an individual track (including daughter cells) which was given a distinct ID number. 
- The links of a given cell at two consecutive time points were listed in links file. 
- The statistics of every cell at every timepoint were listed in spots file, including **x and y coordinates**, **time frame**, and **colour**. These data were analysed by custom-made MATLAB scripts.

## Scritp 1
- To run the script 1, cell tracking files (**links file** and **spots file**) of one experiment were required for input. 

- Users were asked to type in **input file names** as well as an **output file name**. 

- The output file includes values analysed from the cell trajectory, including **track distance**, **track displacement**, **distance/displacement ratio (Tortuosity index)**, **migration velocity** as well as track and step **delta x** and **delta y**. These results were all output in an **Excel file** with several spreadsheets. 

- Noted that the script was run for a single experimental repeat at a time, therefore users would need to change file names accordingly for different repeats.

## Script 2
- Script 2 helped to analyse **cell migration direction** in every step, providing **scatter plots** of cell steps in one control and one tnnt2a repeat, generating a **polar histogram** of cell migration steps (including all control or tnnt2a repeats), as well as the **values** of the polar histogram. 

- This script included two sections.

- The first was to analyse an individual experimental repeat and generate the values of polar histogram in an Excel file. 

- - To run this section, values of **step delta x** and **delta y** were required, which were obtained from **script 1**. 
- - Users were asked to type in **input file names** and **output file names**, **column locations** in the **output** spread sheet for data to be written in and the **number of repeat**. 

- The second section was for plotting, including **scatter plots** and **polar histogram plots**. 
- - User may need to customise the **title** and **colour** of the figures.

## Reference
Tinevez, JY.; Perry, N. & Schindelin, J. et al. (2017), "TrackMate: An open and extensible platform for single-particle tracking.", Methods 115: 80-90, PMID 27713081 (on Google Scholar).

##==============================================================================
##
## Script to run all other scripts
##
## Author: Iman Hosseini-Shakib (ishakib@gmail.com) 
##==============================================================================
## Copyright 2022 Iman Hosseini-Shakib
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================
## Instructions to run:
## This rep is a package of multiple scripts indicated in the order they will be needed. For example, S1_....R indicates step 1.
## The entire package is controlled by main_script.R. This script contains a switch that gives you freedom to run the entire package locally on your machine or use the prepared data that was used in the paper. 
## To use the prepared data, set use_prepared_data=TRUE (this is the default option). To run the code yourself locally, set use_prepared_data=FALSE
## On a regular desktop computer, the entire program if (use_prepared_data=FALSE) should take around 10 hours.
## List of packages that you need to install before running the code is provided below. 
##
## To Run:
##      1. Set the working directory to the main folder (Sobol_Analysis). 
##      2. set use_prepared_data
##      To do so:
##          1. If on RStudio, open the README.md file. Then on the menu bar, go to 
##              Session-->Set Working Directory-->To Source File Location
##          2. If on RStudio, on the lower right box, open "Sobol_Analysis"
##              Then, click on More --> Set as Working Directory
##          3. On console, type the following command: 
##              setwd("~/.../../../../Sobol_Analysis") 
##      3. Run (by clicking on Source in Rstudio)
## What happens next: 
##      R will go through all the scripts one by one   
##      After each script is done, there will be a message om screen reporting that script is done. 
##      The scripts will use the input data saved in the folder called "Inputs"
## Outputs:
##      1. Figures are saved in Figures directory under the main folder
##      2. Data are saved in the Outputs folder under the main directory
## Requirements before running
##      You will need R and these packages: lattice,Kendall,ismev,evdbayes,evir,
##      evd,lhs,fields,plotrix,lhs,rpart,rpart.plot,DEoptim,prim,truncnorm,sdtoolkit,sensitivity,pracma

##==============================================================================
##==============================================================================
##==============================================================================

# Start the program here
rm(list=ls()) #Just in case, remove any variable that is already loaded 




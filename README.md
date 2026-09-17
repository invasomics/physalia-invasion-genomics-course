# Physalia Invasion Genomics Course
## Invasion genomics in (mostly) R

https://www.physalia-courses.org/courses-workshops/invasion-genomics/

### Dr Ang McGaughran
### 12th, 14th and 16th October 2026

### COURSE OVERVIEW
This course on invasion genomics offers a comprehensive three-day program designed to introduce participants to key concepts and techniques in the field. Through a combination of recorded lectures and self-paced tutorials, attendees will explore three main topics: movement dynamics, adaptation dynamics, and genotype-by-environment interactions relevant to invasive species.

On Day1, participants will learn a variety of methods for understanding population diversity and differentiation. Day2 will include an introduction to selection outlier scanning methods. Day3 will focus on the role of environmental variables in shaping adaptive outcomes.

Participants will gain both fundamental knowledge and practical experience, providing a strong starting point for conducting their own research in this rapidly evolving field. They will also have the opportunity to interact with the instructor to discuss questions and address any challenges encountered during the course.

### TARGET AUDIENCE
This course is aimed at graduate students, researchers, and professionals interested in an introduction to invasion genomics. Participants should have a basic understanding of genetics/genomics and evolutionary biology, as well as familiarity with R.

We recommend that participants without prior experience in R complete an introductory R tutorial before attending.

### BACKGROUND READING
• McGaughran A, Rollins LA, Dhami MK (2026). Six recommendations when applying genomic tools to invasion biology. Biological Invasions, 28, 1. https://doi.org/10.1007/s10530-025-03718-x

• McGaughran A, Dhami MK, Parvizi E, Vaughan AL, Gleeson DM, Hodgins KA, Rollins LA, Tepolt CK, Turner KG, Atsawawaranunt K, Battlay P, Congrains C, Crottini A, Dennis TPQ, Lange C, Liu XP, Matheson P, North HL, Popovic I, Rius M, Santure AW, Stuart KC, Tan HZ, Wang C, Wilson J (2024). Genomic tools in biological invasions: Current state and future frontiers. Genome Biology and Evolution, 16, evad230. https://doi.org/10.1093/gbe/evad230

• North HL, McGaughran A, Jiggins C (2021). Insights into invasive species from whole-genome resequencing. Molecular Ecology, 30, 6289–6308. https://doi.org/10.1111/mec.15999

### SESSION CONTENT
During the course days, recorded lectures and self-guided practicals (accompanied by annotated scripts) will be released, allowing participants to progress at their own pace. Live remote support will be provided throughout to accommodate different time zones and assist with any questions.

#### Monday - 12th October (DayOne)
Movement dynamics: inferring invasion sources, dispersal pathways, and population structure.

Recorded lecture: Introduction to invasion genomics, including key research questions, commonly used methods, and analytical challenges associated with invasive species studies.

Tutorial 1: Genetic diversity and differentiation

#### Wednesday - 14th October (DayTwo)
Adaptation dynamics: inferring selection processes in the invaded range.

Recorded lecture: Overview of selection analyses and how they can be used to infer adaptive processes in invasive species and populations.

Tutorial 2: Selection outlier detection methods

#### Thursday - 16th October (DayThree)
Role of the environment: inferring adaptation to local conditions using GEA

Recorded lecture: Introduction to environmental association analysis and its application to understanding adaptation in invasive species.

Tutorial 3: Genotype-by-environment association methods

### COURSE PREPARATION
#### R
The code has been tested on version 4.4.2, so I would suggest that you have at least version 4.2 (and preferably version 4.3 or higher) installed prior to the course start date.
Note that R and RStudio are two different things: updating RStudio will be insufficient, you'll also need to update R.
To check what version of R you have installed, you can run:
```
version
```
The process of updating R will vary depending on your operating system, but you should be able to find instructions for this online.

We will make use of several R packages throughout the course that you'll need to have installed:
- adegenet
- caret
- corrplot
- dplyr
- extendedForest
- geodata
- ggplot2
- ggvenn
- gradientForest
- LEA
- lfmm
- ozmaps
- pcadapt
- qvalue
- raster
- remotes
- reshape2
- robust
- pophelper
- psych
- qvalue
- sp
- StAMPP
- tidyr
- UpSetR
- vegan
- vcfR

You will need to install these packages using a few different methods. Some are Bioconductor packages, which are installed with BiocManager:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("LEA")
BiocManager::install("qvalue")
```
For the package remotes, you'll need to run:
```
library(remotes)
remotes::install_github('royfrancis/pophelper')
```
The rest should work with install.packages in R, for example:
```
install.packages("vcfR")
install.packages("adegenet")
install.packages("StAMPP")
```
Don't underestimate how long it may take to get all of these packages installed! My recommendation is that you install each package one by one, and then confirm the installation has worked by running:
```
library(packagename)
```
For example:
```
library(vcfR)
```
Once you have achieved this, the tutorials will be able to be run. During each tutorial, you will run the library(packagename) command to load the package as required.

#### BayPass
There is one section of the course that is run outside of R, in the program BayPass. So, you will also need to install BayPass on your system as part of the course preparation.
You should follow instructions here for your system:
https://forge.inrae.fr/mathieu.gautier/baypass_public
In my case, I am working on a linux machine and used the following code in a terminal window:
```
git clone https://forge.inrae.fr/mathieu.gautier/baypass_public.git
cd baypass_public/sources
make clean all FC=gfortran
```
You will then need to copy the g_baypass file (or your version of BayPass may be called ifx_baypass or i_baypass) to your working directory (from baypass_public/sources).
You can copy and paste the file across, or do it via the terminal/shell (replacing 'unname' with whatever your's is called):
```
cp /home/uname/baypass_public/sources/g_baypass ~/Desktop/Physalia/DayTwo/
```
For Windows users, here is one option that might work (it involves first installing WSL (Windows subsystem for linux); which enables running a linux environment on your Windows computer without a traditional virtual machine or dual-boot setup), running from the Windows powershell:
```
wsl --install
wsl
wsl --install -d Ubuntu
sudo apt update
sudo apt install build-essential gfortran
cd ~
git clone https://forge.inrae.fr/mathieu.gautier/baypass_public.git
cd baypass_public/sources
make clean all FC=gfortran
cp /home/uname/baypass_public/sources/g_baypass ~/Desktop/Physalia/DayTwo/
```
In case you are unable to get BayPass to run, the relevant output files will be available for you to use during Day Two.

# Usage
The R script cgMLSTcomparison.R aims at performing Principal Component Analysis (PCA) and Generalized Linear Model (GLM) from cgMLST outcomes described in the "Listeria-cgMLST-Additional-file-3.tsv" of the article entitled "In vitro and in silico parameters for precise cgMLST typing of Listeria monocytogenes".
# Dependencies
The R script DownsampledReads.R was prepared and tested with R version 3.6.3 and RStudio 1.3.1093.
- library(ggplot2)
- library(plyr)
- library(reshape2)
- library(tidyr)
- library(AER)
- library(spgs)
- library(usethis)
- library(devtools)
- library(ggbiplot)
# Install R and Rstudio
## 1/ Install R (from configured sources)
```
sudo apt update
sudo apt -y install r-base
R --version
```
## 2/ Install RStudio (from dowloaded rstudio-1.3.1093-amd64.deb)
```
sudo apt install gdebi-core
sudo gdebi /home/Downloads/rstudio-1.3.1093-amd64.deb
rstudio --version
```
# Update R
## 1/ Check current R version
```
R --version
```
## 2/ Update and upgrade apt-get
```
sudo apt-get update
sudo apt-get upgrade
```
## 3/ Check available last recent R version
```
sudo apt-cache showpkg r-base
```
## 4/ Update R
```
sudo apt-get install r-base
```
## 5/ Check updated R version
```
R --version
```
# Install devtools (Ubuntu 20.04)
```
sudo apt install libssl-dev
```
# Follow step per step the R script with RStudio
```
rstudio cgMLSTcomparison.R
```
# Illustration
![PCA figure](https://github.com/Nicolas-Radomski/cgMLSTcomparison/blob/main/illustration.png)
# Acknowledgment
My former colleague Laurent Guillier with whom I learned a lot
# Author
Nicolas Radomski

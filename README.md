<p align="center">

  <h1 align="center">
    [Dowen_Fan]
  </h1>

  <p align="center">
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template" target="_blank">
     <img alt="Status"
          src="https://img.shields.io/badge/status-active-success.svg" />
   </a>
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues" target="_blank">
     <img alt="Github Issues"
          src="https://img.shields.io/github/issues/stjudecloud/bioinformatics-tool-template"  />
   </a>
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/pulls"  target="_blank">
     <img alt="Pull Requests"
          src="https://img.shields.io/github/issues-pr/stjudecloud/bioinformatics-tool-template"  />
   </a>
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/blob/main/LICENSE.md" target="_blank">
     <img alt="License: MIT"
          src="https://img.shields.io/badge/License-MIT-blue.svg" />
   </a>
  </p>


  <p align="center">
   [ChIA-PET data processing and Insulated neighborhood building pipeline] 
   <br />
   <a href="#"><strong>Explore the docs »</strong></a>
   <br />
   <a href="#"><strong>Read the paper »</strong></a>
   <br />
   <br />
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues/new?assignees=&labels=&template=feature_request.md&title=Descriptive%20Title&labels=enhancement">Request Feature</a>
    | 
   <a href="https://github.com/stjudecloud/bioinformatics-tool-template/issues/new?assignees=&labels=&template=bug_report.md&title=Descriptive%20Title&labels=bug">Report Bug</a>
   <br />
    ⭐ Consider starring the repo! ⭐
   <br />
  </p>
</p>

---
## Quick Start
Installing Dowen_Fan:

git clone https://github.com/stjude/Dowen_Fan.git; cd Dowen_Fan

## Description:

Dowen_Fan is a pipeline that enables user to define significantly interacting ChIA-PET loops. This tool is robust and efficient in deducing the Insulated Neighborhoods.

#### Dependencies
cutadapt (v 4.3) 
   <br />
macs14 (v 1.4.2) 
   <br />
bowtie (v 2.3.5.1) 
   <br />
samtools (v 1.15.1) 
   <br />
bedtools (v 2.30.0) 
   <br />
Python 
   <br />
R (v 4.2.2) 
   <br />

## Usage  

sh singlenodeChIA-PETv5.sh -f file_R1.fastq -r file_R2.fastq -b fullpath_bowtieindexname  

where: 

f1 is forward fastq read file from ChIA-PET 

r1 is reverse fastq read file from ChIA-PET 

b is full path including the index name (bowtie indexed genome) 

### Contact
Provide contact information if the author of the code chooses to.  Issues are another effective way of communication as well.

#### Limitations
Are there any know limits to the tool.  Some examples of this might be:
  - Does not run on OSX
  - Version 1.2 of library XYZ causes `this` unexpected behavior
---
#### COPYRIGHT 
Provide St Jude Copyright information here.  Make sure to update the year.

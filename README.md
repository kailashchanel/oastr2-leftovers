# Graph Theoretic Approach to Investigating the Critical Thresholds At Which The Galaxy Filamentary Structure Forms

This is the code for a cosmology study primarily using Julia. This study aims to find the critical density (using percolation theory) that galaxy filaments form at. You can view the paper here: [link to paper, only put when published].

![](https://github.com/AlexanderJCS/JANGL/assets/98898166/0436f41f-7651-492a-8a57-9b154368ba13)
An image of a galaxy cluster in the GAMA datset.


## Data
Some files may reference a file called `GAMA_CZ5Unj.csv`. This is the Data Release 4 of the GAMA (Galaxy and Mass Assembly) survey. You may download the data through [the GAMA website](http://www.gama-survey.org/dr4/).

TODO: provide the specific download link instead of the homepage to GAMA

## Running

This code is structured as a bunch of separate programs. There is no "one program" to run to run the entire analysis. Notes on what each program does can be found in a multi-line comment at the top each file.

TODO: create a bash script that will run all analysis with a single `./run.sh` command

## Requirements

This repository is built with Julia version 1.9.2 in mind.
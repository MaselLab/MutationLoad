# Mutational Drought Project

## Description

This repository is part of the Mutational Drought project and enables the reproduction of results published in Mawass et al. (2024). The project evaluates the role of mutational drought—defined as the shortage of beneficial mutations in shrinking populations—relative to the previously studied mutational meltdown, which involves the accumulation of deleterious mutations.

The repository includes:

- A wholge-genome forward-time evolutionary simulator written in C. To simulate whole genomes efficiently, the genome is divided into discrete linkage blocks. Each block is assigned a numeric value summarizing the effects of all mutations within it. Recombination is limited to hotspots between blocks, with two crossover events per chromosome per meiosis, improving runtime. The simulation is compatible with tskit, enabling detailed tracking and analysis of mutation histories.
- A Mathematica script for analytical analyses and outputting results into CSV dataframes used for plots using the R script.
- Python scripts for analyzing tree sequences output from the simulations.
- R scripts for generating all the plots published in the manuscript.

More details on each component are provided below.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)
- [Contact](#contact)

## Installation
To run the program, first note the following dependencies: 

``` 
tskit -version 1.0.0 (C); gsl – version 2.8 

``` 

# Tskit Download 

To use tskit with the simulation code for tree sequence recording, version 1.0.0 of the tskit GitHub repository should be cloned into your home directory; this location of the repository is assumed by the Makefile for the program. To obtain a path to the home directory, you can use the command  

``` 
echo $HOME 
``` 

on most systems. The following command can also be used move to this home directory: 

``` 
cd $HOME 
``` 

To download the appropriate version of the tskit repository, when in your home directory, the commands given below can be typed into the terminal: 

```
git clone https://github.com/tskit-dev/tskit.git 

cd tskit 

git checkout C_1.0.0
```
# gsl Installation 

To compile and run the simulation, you must install GNU Scientific Library (GSL) version 2.8. You can download it from https://www.gnu.org/software/gsl/, where the package includes an INSTALL file with detailed instructions. GSL should be installed to the default location—do not use the --prefix option to specify a custom directory.

As an example, we provide the installation steps used on one of our systems, though these may vary by operating system. After downloading the .tar file for GSL 2.8, extract it by double-clicking (on some systems) or using the following terminal command:

```
tar –zxvf gsl-2.8.tar.gz 
```
Once the unzipping is done, you should have a new folder called gsl-2.8. Enter the folder with the command 
```
cd gsl-2.8 
```
Next, we configure the installation, which can be done with 
```
./configure 
```
Note that the –prefix option is not used in this case, since we want the folder to be installed in the default location. This command may take a few minutes to complete, and once it’s done, the library can be compiled with 
```
make 
```
Again, this step will require a few minutes to complete. To ensure that the compilation was completed properly, we run the following command, checking to make sure that no errors arise: 
```
make check 
```
Finally, once this command runs, the program can be installed with 
```
make install 
```
With this, gsl should have been installed to the default location in your directory system!

# Compilation 
Once all of the necessary dependencies have been downloaded and installed, the code can be compiled using the Makefile within the repository. First, though, the HOME_DIR macro at the top of the Makefile should be adjusted so that it contains the name of your own home directory, which will enable make to find the files necessary for compilation. Again, to get the path to the home directory, the command given below works on most systems: 

```
echo $HOME 
``` 
Running the make command will then compile the program, creating an executable named mutationload that can be used to run the simulation itself:
```
make
```

## Usage

- The simulator can be compiled using the Makefile or executed via provided Bash scripts in an HPC Slurm environment.

# Running the Program 

Considering that the simulation takes in many command-line arguments, a bash shell script is used to run the executable; more specifically, the bash file titled “bash_local.sh” within the repository acts as a good starting place for such a shell script. Make sure, though, to first use the chmod command to make both the mutationload file and the bash script executable:
```
chmod +x bash_local.sh
chmod +x mutationload
```
Once this has been done, we can invoke our bash script to run the mutationload executable:

``` 
./bash_local.sh 
``` 

The script can then be adjusted to change the parameters for a particular simulation, with the meanings of the bash parameters being described as follows: 

### Relative Fitness Parameters 
+ timeSteps – number of generations the simulation should run for 

+ initialPopsize –initial number of individuals in the population 

+ mud – genome-wide deleterious mutation rate (Ud) 

+ chromosomesize – number of linkage blocks per chromosome 

+ numberofchromosomes – number of distinct chromosomes per individual (note that individuals are diploid, so they have two copies of each chromosome) 

+ bentodelratio -- ratio of Ub (beneficial mutation rate) to Ud (deleterious mutation rate) 

+ sb – mean effect of a beneficial mutation

+ bendist – indicates the distribution for beneficial effect size; 0 for point, 1 for exponential, and 2 for uniform 

+ typeofrun – indicates the type of simulation run; 0 for root for Sb, 1 for single, 2 for root for Ncrit

+ slope – indicates the slope value of the fitness flux to target when finding root

+ seed – seed value for RNG

+ K – population capacity (specific for absolute model)

+ fitnesstype – 0 for a run using relative fitness, 1 for a run using absolute fitness 

+ r – strength of global epistasis (usually set for 1 to exclude epistasis)

+ i_init – any value (specific for absolute model)

+ s – any value (specific for absolute model)

+ tskitstatus – 0 for no tskit, 1 to have tskit on the entire simulation, 2 to have tskit on only after the burn-in phase

+ modularepis – 0 for no modularity, 1 for modular epistasis

+ elementsperl – value of modular elements per linkage block (only valid if modularepis=1)

+ snapshot – 0 to not use snapshot of previous run, or 1 to use snapshot or previous run

+ file1 – snapshot file

+ SdtoSbratio – ratio of Sd (deleterious effect size) to Sb (beneficial effect size). If deldist=0, this parameter is neglected. 

+ deldist – determines the type of distribution for deleterious effect size; 0 to use the distribution inferred for humans in Kim et al. 2017, 1 for exponential, 2 for point

+ rawdatafilesize – controls how many lines to print out in the raw data file (helps compact the size of the file if timeSteps is large)

+ redinmaxpopsize – 0 to exclude reduction in max pop size, 1 to include reduction in pop size (specific for absolute model)

+ calcfixation – 0 to not perform fixation time calculation using trees, 1 to calculate fixation time (specific for absolute model; NOTE: incomplete feature)

- The Mathematica script requires **Mathematica Desktop** or **Mathematica Online**.
- The Python script can be executed in a bash terminal in a conda environment(see script headers for dependency details).

Use the following line to activate the conda environment for the required dependencies to run the python script
```
cd $HOME
conda env create -f tskit_env.yml
conda activate tskit_env
```

- The R scripts can be run via Bash or **RStudio** (load R environment using renv file for dependencies).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For inquiries, please contact:

- **Email:** [masel@arizona.edu](mailto:masel@arizona.edu)

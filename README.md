 # Simulation Overview 

For the purposes of studying topics in theoretical population genetics, this repository contains a whole-genome forwards time simulation constructed by the Masel lab at the University of Arizona. To simulate whole genomes efficiently, this simulation breaks the genome up into discrete linkage blocks, with each of these being associated with a numeric value summarizing the combined effects of all mutations occurring within that linkage block. Recombination then only occurs at hotspots between adjacent linkage blocks, with two recombination events occurring per chromosome per meiosis, allowing a reasonable runtime for the simulation. Our code is then compatible with tskit, which allows us to track and analyze the history of individual mutations within the population. 

 

# Dependencies 

To run the program, first note the following dependencies: 

``` 

tskit – version 0.5.8; gsl – version 2.8 

``` 

# Tskit Download 

To use tskit with the simulation code for tree sequence recording, version 0.5.8 of the tskit GitHub repository should be cloned into your home directory; this location of the repository is assumed by the Makefile for the program. To obtain a path to the home directory, you can use the command  

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

	git checkout beafeba3f576da4524e24b00aa184e608f6de4c4 

``` 

 

# gsl Installation 

Version 2.8 of the GNU Scientific Library, or gsl, must also be installed to compile and run the simulation code. Information on this library and how to download it on your own system can be found here (https://www.gnu.org/software/gsl/), and the downloaded package will then contain an INSTALL text file, providing specific information on how to install gsl. For the purposes of running the simulation, gsl should be installed to the default location, rather than using the --prefix option to specify a particular folder.  

 

To get an idea of how this installation process will work, we provide the set of steps used to install this library on one of our own systems, although this exact approach likely won’t work for all operating systems. Downloading version 2.8 of gsl from the site linked above should provide you with a .tar file, which must then be unzipped. On a local system, this can often be done by double-clicking the .tar within your file navigator; otherwise, for most systems, the following tar command should work: 
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
Note, though, that if using the older version of the Makefile for this repository, titled "old_makefile", compilation can be performed by specifyng this alternate file:
```
	make -f old_makefile
```
Also, if an error comes up during compilation due to a missing file, this is frequently caused by an issue related to the location of the "dependencies" folder within the simulation's directory. It can be helpful to make sure that the location of this folder matches that which is specified in the Makefile (i.e. that the folder is within the larger MutationLoad folder). If issues persist, it's also possible to try moving all of the files within the folder into the overall MutationLoad folder, then removing mentions of the "dependencies" folder from the Makefile. 
 

# Running the Program 

Considering that the simulation takes in many command-line arguments, a bash shell script is used to run the executable; more specifically, the bash file titled “general_bash_local.sh” within the repository acts as a good starting place for such a shell script. Make sure, though, to first use the chmod command to make both the mutationload file and the bash script executable:
```
	chmod +x general_bash_local.sh
	chmod +x mutationload
```
Once this has been done, we can invoke our bash script to run the mutationload executable:

``` 

	./general_bash_local.sh 

``` 

The script can then be adjusted to change the parameters for a particular simulation, with the meanings of the bash parameters being described as follows: 

### Relative Fitness Parameters 

+ fitnesstype – 0 for a run using relative fitness, 1 for a run using absolute fitness 

+ timeSteps – number of generations the simulation should run for 

+ initialPopsize –initial number of individuals in the population 

+ mud – genome-wide deleterious mutation rate (Ud) 

+ chromosomesize – number of linkage blocks per chromosome 

+ numberofchromosomes – number of distinct chromosomes per individual (note that individuals are diploid, so they have two copies of each chromosome) 

+ bentodelratio -- ratio of Ub (beneficial mutation rate) to Ud (deleterious mutation rate) 

+ sb – mean effect of a beneficial mutation

+ bendist – indicates the distribution for beneficial effect size; 0 for point, 1 for exponential, and 2 for uniform 

+ typeofrun – indicates the type of simulation run; 0 for root, 1 for single 

+ tskitstatus – 0 for no tskit, 1 to have tskit on the entire simulation, 2 to have tskit on only after the burn-in phase  

+ deldist – determines the type of distribution for deleterious effect size; 0 to use the distribution inferred for humans in  Kim et al. 2017, 1 for exponential, 2 for point 

+ SdtoSbratio – ratio of Sd (deleterious effect size) to Sb (beneficial effect size). If deldist=0, this parameter is neglected. 

 

### Absolute Fitness Parameters 

+ slope – the slope for the contour line 

+ seed – the seed for the random number generators in the program 

+ K – the carrying capacity in a deterministic population with minimal load 

+ r – epistasis variable R, which varies between 0 and 1, with 1 representing no epistasis and 0 representing full epistasis 

+ i_init – the number of mutations away from the worst-case scenario at the beginning of the simulation 

+ s – selection coefficient of a beneficial mutation at minimum death rate (i.e. max birth rate) 

+ rawdatafilesize – raw data file size in data points, which determines the sampling rate of data points 

+ redinmaxpopsize – input for reduction in max population size 

+ calcfixation – 0 to have fixation calculation off, 1 to have fixation calculation on  

 

### Modular Epistasis Parameters 

+ modularepis – 0 for runs without modular epistasis, 1 for runs with modular epistasis. Note that modularepis=1 works only in absolute fitness runs (fitnesstype=1) 

+ elementsperl – number of elements per linkage block 

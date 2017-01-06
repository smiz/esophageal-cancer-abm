# esophageal-cancer-abm

See main.cpp for details about the model.

To run this model, take the following steps:

(1) Create an input data file.

An example data file is input.txt. The parameters are
 diffusion_rate    - The rate of spread in mm^2/year of cancer and dysplasia
 be_onset_age      - The mean age in years at which BE appears
 mutate_be         - The number of mutations per cell per year into dysplasia
 mutate_dysplasia  - The number of mutations per cell per year into cancer
 stem_cell_density - The number of stems cells per mm^2 in the BE segment

(2) Run the program with its command line arguments. Here are some examples.

Run the simulation using random seed 10, the input file input.txt, and 
create a snapshot for biopsy at ages 30, 40, and 50 years.

 ./a.out -ranseed 10 30 40 50

Run the simulation using the default random seed, the input file mine.txt, and 
create a snapshot for biopsy at ages 40 and 50 years.

 ./a.out -var mine.txt 40 50

Run the simulation using random seed 2, the input file mine.txt, and 
create a snapshot for biopsy at ages 40 and 50 years.

 ./a.out -ranseed 2 -var mine.txt 40 50

(3) Look at the output.

At each biopsy instant, a count of cell types will be printed to the screen.
A file tumor.csv.0 will be created for the first biopsy, tumor.csv.1 for the
second, etc. The format of these files is a header followed by a list of
locations that have cancer or dysplasia. Here is an example of a file with
cancer at point 10, 10, 10 and dysplasia at 11, 10 10

xcoord,ycoord,zcoord,type
10,10,10,3
11,10,10,2

These files can be visualized using paraview.

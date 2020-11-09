# CrisprCountsAnalysis (CCA)

Please take a look at the documentation in the doc subdirectory.   They will
contain the input file format and more explanations for running CCA.


-------------------------------------------------------------------------------


To build the image where you downloaded the source, do

docker build -t cca .

in the directory with the Dockerfile.  This is optional and one can instead
pull the Docker image from the Docker hub by doing

docker pull tohsumirepare/cca


-------------------------------------------------------------------------------


To run the program, please read the documentation in the docs directory.  Then,
once example files EX1.txt and EX1.repmap are ready in directory, for example, 
/Users/login/temp, then do

docker run --rm -it -v "/Users/login/temp:/home:rw" tohsumirepare/cca crisprCountsAnalysis COUNTS=EX1.txt REPMAP=EX1.repmap OUT_MAX=EX1

in the directory with the EX1.* files if you have pulled the Docker image.  If 
you have compiled the code instead then do


docker run --rm -it -v "/Users/login/temp:/home:rw" cca crisprCountsAnalysis COUNTS=EX1.txt REPMAP=EX1.repmap OUT_MAX=EX1




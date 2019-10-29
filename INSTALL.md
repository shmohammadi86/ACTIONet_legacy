# Installation guideline

## Installing igraph c++ library (needed for running leiden algorithm)
If you encounter the following error during the installation, it means that igraph C++ library is not installed on your system (yet!):

*../inst/include/leiden/GraphHelper.h:4:27: fatal error: igraph/igraph.h: No such file or directory*

Here is how you install it.

### Ubuntu


### Compile from scratch
If you do not have root access, you need to compile igrpah from scratch, then install ACTIONet manually:

** Compiling igraph **

```{bash}
# Compiling igraph from source
mkdir workspace && cd workspace
mkdir -p libs/igraph
wget https://igraph.org/nightly/get/c/igraph-0.7.1.tar.gz
tar -xzvf igraph-0.7.1.tar.gz
cd igraph-0.7.1
./configure --prefix=`pwd`/../libs/igraph/ && make -j8 && make install

# Installing ACTIONet manually
cd ..
wget https://github.com/shmohammadi86/ACTIONet/archive/R.zip 
unzip R.zip
export LD_LIBRARY_PATH=`pwd`/libs/igraph/lib
cd ACTIONet-R/
echo "PKG_CXXFLAGS+=-I../../libs/igraph/include/" >> src/Makevars
echo "PKG_LIBS+=-L../../libs/igraph/lib -Wl,-rpath,`pwd`/../../libs/igraph/lib" >> src/Makevars
R CMD INSTALL .
```





### Linux

Install (upgrade to) the latest version of gcc compiler:

```{bash}
sudo add-apt-repository ppa:jonathonf/gcc-9.0
sudo apt-get install g++-9
```



Install igraph c++ interface

```{bash}
sudo apt-get install libigraph0-dev
```

Or install it from the source from https://igraph.org/c/.



Create  *~/.R* folder and copy *Makevars_gcc.txt* there under the name *Makevars*: 

```{bash}
 mkdir -p ~/.R
 wget --no-check-certificate https://raw.githubusercontent.com/shmohammadi86/ACTIONet/master/Makevars_gcc -O ~/.R/Makevars
```

Install Bioconductor dependencies:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install('scran', 'scater')
```

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```



### Mac OS

Install (upgrade to) the latest version of gcc compiler:

```{bash}
brew install gcc@9
```

Install igraph c++ interface

```{bash}
brew install igraph
```

Or install it from the source from https://igraph.org/c/.

Create  *~/.R* folder and copy *Makevars_gcc.txt* there under the name *Makevars*: 

```{bash}
 mkdir -p ~/.R
 wget --no-check-certificate https://raw.githubusercontent.com/shmohammadi86/ACTIONet/master/Makevars_gcc -O ~/.R/Makevars
```

Install Bioconductor dependencies:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install('scran', 'scater')
```

Now you can install ACTIONet package with devtools:

```{r}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet", ref = "R-release")
```


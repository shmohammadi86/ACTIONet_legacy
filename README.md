# Main ACTIONet repository

# Installation 

## Linux

Install igraph c++ interface -- in terminal

```{bash, echo = TRUE}
sudo apt-get install libigraph0-dev
#Or install it from the source from https://igraph.org/c/.
```

Now you can install ACTIONet package with devtools in R:

```{r, echo = TRUE}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet")
```

## Mac OS

Install the latest version of clang compiler with OpenMP support -- in terminal

```{bash, echo=TRUE}
brew install llvm
```

Install igraph c++ interface -- in terminal

```{bash, echo=TRUE}
brew install igraph
#Or install it from the source from https://igraph.org/c/.
```

Make sure updated clang is used for compiling R packages -- in terminal

```{bash, echo=TRUE}
mkdir -p ~/.R && echo -e "CXXFLAGS=-std=c++14 -lstdc++ -w -m64 -fPIC -march=native -O4\nCXX=/usr/local/opt/llvm/bin/clang++\n" > ~/.R/Makevars
```

Now you can install ACTIONet package with devtools in R:

```{r, echo = TRUE}
install.packages("devtools")
devtools::install_github("shmohammadi86/ACTIONet")
```



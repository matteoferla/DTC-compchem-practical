# Installation

## Overview
### Drives
First and foremost, be aware that there are two drives in your system, one is the local disc and the other is your mounted networked drive. All operations should be in the former to avoid latency.

### Conda
We will be using `conda` for our Python. This is a tool that allows one to make a playground where we can install things without being root and if we break it we can start again! This has these three key subcommands:

* `conda install`
* `conda create` — create a virtual environment
* `conda activate` / `conda deactivate` — turn on a virtual environment


### Jupyter notebooks
We also will be using Python in a Jupyter notebook. A Jupyter notebook is a web-based way to interact with a Python (or Julia) kernel. Namely, you run the webapp serving on a given port, which you visit in your browser.
A variant of this is [Colab](https://colab.research.google.com/), which is an offshoot of an early Jupyter notebook and served off Google's cloud infrastructure.
Vanilla Jupyter has evolved into jupyterlab, which has better features and will be used in this teaching lab.
For more info about Jupyter see [Jupyter notes](jupyter.md)

## Installing conda
Let's say you have the path of the former as `$WORKINGDIR_PREFIX`:

```bash
WORKINGDIR_PREFIX='/👾👾👾/👾👾👾'
```

Open the terminal and install conda.
Change the `$CONDA_ROOT` variable to whatever you want, here the hidden folder `.conda` to be tidy,
but can be whatever you like.

If not already present then skip this block.

Footnote on environment variables: Shell variables are declared without the dollar sign and with no space between the name and the equals sign.
To get a variable use a dollar sign and make sure the letters do not merge w/ the next word,
to avoid this use curly brackets, e.g. `${FOO}1`. If you want to run a command with given environment variable,
then prefix the command like this `FOO=1 BAR=2 bash test.sh`.

```bash
CONDA_ROOT=$WORKINGDIR_PREFIX/.conda
# ---------- download ---------- 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# ----------  install ----------
# The following installs in location specified w/ argument `-p`
# The argument `-b` says to install it without asking questions.
# This does not automatically enable conda, hence the next few lines.
# If you want to always activate conda, add  the lines below this `turn conda on` to your`.bashrc`.
bash Miniconda3-latest-Linux-x86_64.sh -p $CONDA_ROOT -b
# ---------- enable conda ---------- 
__conda_setup="$($CONDA_ROOT'/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f $CONDA_ROOT"/etc/profile.d/conda.sh" ]; then
        . $CONDA_ROOT"/etc/profile.d/conda.sh"
    else
        export PATH=$CONDA_ROOT"/bin:$PATH"
    fi
fi
unset __conda_setup
unset CONDA_ROOT
```
If you have installed conda correctly your prompt should say "(base) $".
Base is not a virtual environment.
If it does not and you just installed it without `-b` argument in the bash then
`source .bashrc` or `.zshrc` if your OS is a MacOS (i.e. you aren't using a teaching lab machine).

The command `conda` will now be callable and the following environment variable present:

* `$CONDA_PREFIX`
* `$CONDA_EXE`

 '': '/data/xchem-fragalysis/mferla/conda/bin/conda',
 '_CE_CONDA': '',
 'CONDA_PREFIX_1': '/data/xchem-fragalysis/mferla/conda',
 '': '/data/xchem-fragalysis/mferla/conda/envs/biochem38',
 'CONDA_PYTHON_EXE': '/data/xchem-fragalysis/mferla/conda/bin/python',
 'CONDA_DEFAULT_ENV': 'biochem38'}

## Creating a virtual environment w/ Jupyter installed
Now let's make a conda virtual environment.

```bash
# ---------- Install goodies ----------
# nb_conda_kernels makes all conda kernels usable in a notebook
# mamba is a fix to conda for speed
# -y stops it asking Qs
# -c is a channel
#  -n base  tells it to install to the base
conda update -n base -y -c defaults conda
conda install -n base -y -c conda-forge nb_conda_kernels mamba
# ---------- make conda env ----------
# `conda create` or `mamba create` creates a virtual environment
# often you have a non-updated OS and the GNU library for C isnt the latest, 
# hence the `CONDA_OVERRIDE_GLIBC` variable below
# this will a  python=3.9 environment
CONDA_OVERRIDE_GLIBC=2.36 mamba create -n compchem -y -c conda-forge python=3.9 nodejs jupyterlab
```
To avoid copy pasting tokens for permission let's set a password (warning this is not hashed):
```bash
jupyter notebook password
```
Now if you run `jupyter-lab` your browser will pop up on page `localhost:8888/` asking for a password.
But first lest install some more things.

## PyRosetta
We will be using a package that require PyRosetta.
This requires registration which will provide with universal username and password.
Matteo wrote a helper module to install this:

```bash
yes | pip install pyrosetta_help
PYROSETTA_USERNAME=👾👾👾 PYROSETTA_PASSWORD=👾👾👾 install_pyrosetta
```

However, for speed within the Stats department (behind `statistics.frodo.ox.ac.uk` (`172.24.83.17`)) today you can instead do:

```bash
yes | pip install https://www.stats.ox.ac.uk/~ferla/pyrosetta-2022.46+release.f0c6fca0e2f-cp39-cp39-linux_x86_64.whl
```

## Other packages
```bash
conda install -c conda-forge trident-chemwidgets
jupyter labextension install trident-chemwidgets
yes | pip install Fragmenstein gist-import rdkit-pypi rdkit-to-params biopython pyrosetta-help michelanglo-api nglview smallworld-api pandas plotly

conda clean --tarballs -y
conda clean --source-cache -y
pip cache purge
```

yes pipe is to not press the Any key to say yes (or whatever word you add after it)
—bonus: in a Dixons shop open the terminal on a demo Mac and troll the shop assistant...

Note that if you are using Colab (ie. not following instructions) a different widget using JSME is required:
```
!pip install git+https://github.com/matteoferla/JSME_notebook_hack.git
```

## Copy the notebooks in this repo

```bash
git clone https://github.com/matteoferla/DTC-compchem-practical.git
```

## Start Jupyterlab

Go to the folder you would want to be root and run:

```bash
jupyter-lab
```

Go to [localhost:8888](http://localhost:8888) and open to `1-redocking.ipynb`.
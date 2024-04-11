# De novo motif discovery

## Installation

This software requires Python 3.0 <=.
It also requires the ability to install packages, with the likes of pip.

To download the files in this repository, please run:

```
git clone https://github.com/nicolayro/motif-discovery
cd motif-discovery
```

Numpy is required for the main program. It can be installed with pip:

```
pip install numpy
```


Creating the sequence logo requires Biopython as well, but it is not part of the main script:
```
pip install biopython
```


## Usage

To run the tool, use the following command:

```
python main.py <filename> <motif_width>
```

An example with the supplied test data set:

```
python main.py crp0.fna 18
```

The script can also be run without arguments. In this case a very small test data set is used.

```
python main.py
```

## Visualizing the sequence logo

The main script create an `output.jaspar` file. If you wish to visualize this data, you can create a sequence logo with the command:

```
python visualize.py output.jaspar
```


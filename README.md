# PhaKinPro

PHArmacoKINetic PROperty calculator: Used to predict varies pharmacokinetic properties using QSAR models. If you use please cite [our paper](https://pubs.acs.org/doi/10.1021/acs.jmedchem.3c02446). There is a [webserver](http://34.170.18.221/) that runs these models, but for large numbers of compounds, running locally using this code is much more effective

## Installation guide

### Prerequisites

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html)

### Create environment

```bash
git clone https://github.com/jonswain/PhaKinPro.git
cd PhaKinPro
conda env create -f environment.yml
activate PhaKinPro
```

## Example usage

The packages necessary to run the project are now installed inside the conda environment.

After downloading, `PhaKinPro/phakinpro.py` can be called from the command line with `python phakinpro --help`

`--infile` is required and is the fileloc for a csv of SMILES to predict properties for. Requires that csv has header and is comma seperated
`--smiles_col` is the name of the column containing the SMILES strings of interest. Defaults to "SMILES"
`--id_col` is the name of the column containing the compound identifiers. Defaults to "identifer"
`--outfile` is the fileloc of where the output csv file should go. Defaults to `current-working-directory/phakin_output.csv`
`--ad` flag to turn on applicability domain calculation for the models

```bash
python PhaKinPro/phakinpro.py --infile input.csv --smiles_col smiles --id_col indentifiers --outfile out.csv
```

## Webserver interface

This repository also contains the code to run a local webserver (or host your own). You can start the server by running `qunicorn wsqi:app` (or using the devolpment flask server by setting the `FLASK_APP` variable: `$env:FLASK_APP = "main"` on windows or `export FLASK_ENV=main` on unix). From that access 127.0.0.1:5000 to view the local server

Thanks to JSME for a free and easy to use molecule editor for webpages
[Bienfait, B., Ertl, P. JSME: a free molecule editor in JavaScript. J Cheminform 5, 24 (2013).](https://doi.org/10.1186/1758-2946-5-24)

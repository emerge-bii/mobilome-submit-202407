# Setup

## Create conda environment

```bash
git clone https://github.com/emerge-bii/mobilome-submit-202407.git
cd mobilome-submit-202407
mamba env create -f requirements.yaml --force
conda activate mobilome_submit_v0*
```

## download data package

- Option 1: You can download the data package from this url (under contructon) and
name it as `som-data`.
- Option 2: If you are on OSC, it's already there under the EMERGE globus directory. You can make symlink named `som-data` by `ln -s /fs/ess/PAS1212/EMERGE_sharing/mobilome_paper/som-data`

# Reproduce figs

```bash
./generate-fig.sh
```

You can see figs main figures: Fig.1, Fig.2, Fig.4 and supplementary figures: fig.S2 - fig.S12 in `fig.outdir`. Note that Fig.3, Fig.5, fig.S1 were manually made so not included here.

Test data for umccrise
----------------------------------------------------------------------------

## Generating test data

Install dependencies for generating test data:

```
conda install -c bioconda pybedtools sambamba tabix snakemake-minimal
```

Run on a bcbio-nextgen WGS project:

```
snakemake -p Snakefile.prep_test_data --config final=/data/cephfs/punim0010/data/Results/Patients/CUP_SC932/final out=data/bcbio_test_project
```

Or on RNAseq WTS project:

```
snakemake -p Snakefile.prep_test_data_rnaseq --config final=final out=data/bcbio_test_project_rnaseq
```

## Testing

```
nosetest -s test.py
```
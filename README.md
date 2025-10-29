# plyCheck

Annotate ply alleles

``db/ply_alleles.txt``: allele sequences extracted from NCBI

``db/ply_alleles.fasta``: allele sequences to be used in the pipeline

Currently known alleles are [20](https://doi.org/10.1371/journal.pone.0134055), but the sequences for ply-14 and ply-20 are not available yet in the db.

# Usage

Prepare input files in ``config/config.yml`` file.

Create output directories ``òut/blast`` and ``out/logs``

Then  run snakemake as:

``snakemake --use-conda --cores 4 -R all -p``

# Outputs

| output  | description |
|----------|----------|
| ``blastoutput.tsv``    | Merged output from blast with statistics and initial ply allele assignation    |
| ``alleles.fna``    | File with extracted ply genes    |
| ``aa_unaln.fasta``    | File with extracted ply aa and reference aa   |
| ``aa_aln.fasta``    | Aligned  ``aa_unaln.fasta``    |
| ``aa_aln_corrected.fasta``    | Corrected ply allele assignment    |
| ``aa_changes.tsv``    | Initial ply allele assignment    |
| ``aa_changes_corrected.tsv``    | Corrected ply allele assignment    |

To generate the file below:
``python workflow/scripts/update_fasta_headers.py --input-fasta out/aa_aln.fasta --mapping-file out/aa_changes_corrected.tsv --output-fasta out/aa_aln_corrected.fasta``
The ``aa_aln_corrected.fasta`` can be used with IQTREE to generate a phylogenetic tree as in:

``iqtree -s out/aa_aln_corrected.fasta -m MFP -bb 1000``

# Author
Judit Burgaya, judit.burgaya@gmail.com | judit.burgayaventura@unibe.ch

## Table of contents

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Extract the CheckM stats for an assembly accession](#extract-the-checkm-stats-for-an-assembly-accession)
- [Get collection date and geographical location for a list of biosample accessions](#get-collection-date-and-geographical-location-for-a-list-of-biosample-accessions)
- [Convert a VCF file to BED](#convert-a-vcf-file-to-bed)
- [Download an assembly, or assemblies, by accession](#download-an-assembly-or-assemblies-by-accession)
- [Change the chromosome name in a VCF](#change-the-chromosome-name-in-a-vcf)
- [Get depth of coverage as a BED from a BAM file](#get-depth-of-coverage-as-a-bed-from-a-bam-file)
- [Sort a fastq file by length without reading the whole thing into memory.](#sort-a-fastq-file-by-length-without-reading-the-whole-thing-into-memory)
- [Download summaries for all bacterial assemblies in RefSeq](#download-summaries-for-all-bacterial-assemblies-in-refseq)
- [Get Run accessions for a BioSample accession](#get-run-accessions-for-a-biosample-accession)
- [Install and use Aspera to download from ENA/SRA](#install-and-use-aspera-to-download-from-enasra)
- [Extract taxon ID for list of GenBank accessions](#extract-taxon-id-for-list-of-genbank-accessions)

<!-- TOC end -->


---

<!-- TOC --><a name="extract-the-checkm-stats-for-an-assembly-accession"></a>
### Extract the CheckM stats for an assembly accession

```
$ datasets summary genome accession GCF_019454365.1 | jq '.reports' | jq '.[0]' | jq '.checkm_info'
{
  "checkm_marker_set": "Mycobacterium tuberculosis",
  "checkm_marker_set_rank": "species",
  "checkm_species_tax_id": 1773,
  "checkm_version": "v1.2.2",
  "completeness": 95.36,
  "completeness_percentile": 1.6680924,
  "contamination": 0.37
}
```

or from python

```python
import requests
import json

acc = "GCF_019454365.1"
url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/dataset_report"
response = requests.get(url)
if response.status_code == 200:
    data = response.json()
    report = data["reports"][0]
    checkm_info = report["checkm_info"]
else:
    print("Error:", response.status_code)
    
print(json.dumps(report["checkm_info"], indent=4))
{
    "checkm_marker_set": "Mycobacterium tuberculosis",
    "checkm_species_tax_id": 1773,
    "checkm_marker_set_rank": "species",
    "checkm_version": "v1.2.2",
    "completeness": 95.36,
    "contamination": 0.37,
    "completeness_percentile": 1.6680924
}
```

---

<!-- TOC --><a name="get-collection-date-and-geographical-location-for-a-list-of-biosample-accessions"></a>
### Get collection date and geographical location for a list of biosample accessions

```
efetch -db biosample -input accessions.txt -mode xml | 
xtract -pattern BioSample -element BioSample@accession \
    -block Attribute -if Attribute@attribute_name -equals "collection_date" \
    -or Attribute@attribute_name -equals "geo_loc_name" \
    -sep ": " -element Attribute@attribute_name,Attribute
```

which outputs 

```
SAMN35995983    collection_date: 2022   geo_loc_name: USA
SAMN35995984    collection_date: 2017   geo_loc_name: USA
```

---

<!-- TOC --><a name="convert-a-vcf-file-to-bed"></a>
### Convert a VCF file to BED

```
bcftools query -f '%CHROM\t%POS0\t%POS\n' in.vcf
```

---

<!-- TOC --><a name="download-an-assembly-or-assemblies-by-accession"></a>
### Download an assembly, or assemblies, by accession

```
ncbi-genome-download -A GCF_000285655.3 -P -F fasta bacteria
```

`-A` can be a comma separated list of a file to accessions

---

<!-- TOC --><a name="change-the-chromosome-name-in-a-vcf"></a>
### Change the chromosome name in a VCF

echo -e "old_chr\tnew_chrom" | bcftools annotate --rename-chrs - in.vcf

---

<!-- TOC --><a name="get-depth-of-coverage-as-a-bed-from-a-bam-file"></a>
### Get depth of coverage as a BED from a BAM file

```
$ bedtools genomecov -ibam in.bam -bga > out.bed
```

`-bga` provides the output as a [bedgraph](https://genome.ucsc.edu/goldenPath/help/bedgraph.html) file that also reports zero-depth sites.

The output will look like this

```
NC_000962.3     0       35      29
NC_000962.3     35      80      30
NC_000962.3     80      106     31
NC_000962.3     106     107     30
NC_000962.3     107     126     31
```

For per-base depth, prefer

```
$ bedtools genomecov -ibam in.bam -d
```

or 

```
$ mosdepth prefix in.bam
```

---

<!-- TOC --><a name="sort-a-fastq-file-by-length-without-reading-the-whole-thing-into-memory"></a>
### Sort a fastq file by length without reading the whole thing into memory.

```
gzip -dc in.fq.gz | 
  paste - - - - | 
  perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;' | 
  sort -n | 
  cut -f2- | 
  tr "\t" "\n" > len_sorted.fq
```

Credit: <https://thegenomefactory.blogspot.com/2012/11/sorting-fastq-files-by-sequence-length.html>.

This works surprisingly fast and with very little memory. On a 431 MB (compressed) fastq file, it took 7.27 seconds (wall clock) and used 1.2 MB of memory.

You can then sanity check the output to see if it is sorted in descending (`-r`) or ascending (default) order with this python script

```python
"""This script checks if a fastq file is sorted by length"""

import sys
import os
import argparse


def is_sorted(fastq_file, reverse=False):
    """Check if a fastq file is sorted by length"""

    with open(fastq_file) as f:
        last_len = 0 if not reverse else float("inf")
        for i, line in enumerate(f, start=1):
            if i % 4 == 2:
                seq_len = len(line.strip())
                if reverse:
                    if seq_len > last_len:
                        return False
                else:
                    if seq_len < last_len:
                        return False
                last_len = seq_len

    return True


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Check if a fastq file is sorted by length"
    )
    parser.add_argument("fastq_file", help="fastq file to check", default=0, nargs="?")
    parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        help="check if the file is sorted in descending (reverse) order",
    )
    args = parser.parse_args()

    if args.fastq_file not in [0, "-"] and not os.path.exists(args.fastq_file):
        print("Error: file not found")
        sys.exit(1)

    if args.fastq_file == "-":
        args.fastq_file = 0

    if is_sorted(args.fastq_file, reverse=args.reverse):
        print("The file is sorted by length")
        sys.exit(0)
    else:
        print("The file is NOT sorted by length")
        sys.exit(1)


if __name__ == "__main__":
    main()

```

Which can be used as 

```
$ python is_len_sorted.py len_sorted.fq
The file is sorted by length
```
---

<!-- TOC --><a name="download-summaries-for-all-bacterial-assemblies-in-refseq"></a>
### Download summaries for all bacterial assemblies in RefSeq

Using the NCBI's [`datasets`](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) command line tools.

```
$ datasets summary genome taxon bacteria --as-json-lines --assembly-source RefSeq --assembly-level complete --assembly-version latest --mag exclude > bacteria.jsonl
```

this ensures that we only get summaries for 'complete' assemblies

> A Complete assembly in RefSeq means that the genome is fully assembled, with no gaps, and typically represents the entire genome. All chromosomes and extra-chromosomal elements (e.g., plasmids or organelle genomes) are completely represented.

This can then be converted to a TSV and extract fields of interest with

```
FIELDS="accession,organism-name,organism-tax-id,assminfo-bioproject,assminfo-biosample-accession,assminfo-name"
dataformat tsv genome --inputfile bacteria.jsonl --fields $FIELDS > bacteria.tsv
```

a full list of fields can be found [here](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/tsv/dataformat_tsv_genome/).

You can also download these summaries for a single accession with

```
datasets summary genome accession GCF_000006945.2 | jq '.reports[0].assembly_info.biosample.accession'
```

the above example will give you the biosample accession.

---

<!-- TOC --><a name="get-run-accessions-for-a-biosample-accession"></a>
### Get Run accessions for a BioSample accession

```
$ biosample=SAMN31564381
$ fields=run_accession,instrument_platform
$ curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&query=sample_accession=${biosample}&fields=${fields}"
run_accession   instrument_platform
SRR22225500     ILLUMINA
SRR22225499     OXFORD_NANOPORE
```

the full list of available fields can be found [here](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/advanced-search.html#what-fields-can-i-use-in-my-search). Alternatively, the [EBI advanced search](https://www.ebi.ac.uk/ena/browser/advanced-search) can be used to construct queries and then copy the corresponding `curl` command - and there is also a link to the API docs. The swagger API docs where you can run test queries is [here](https://www.ebi.ac.uk/ena/portal/api/swagger-ui/index.html).

Alternatively, if there are no runs, you can use the sample data type as such

```
curl "https://www.ebi.ac.uk/ena/portal/api/search?result=sample&format=tsv&query=sample_accession=${biosample}&fields=all"
```

To do the inverse, that is search for a BioSample accession from a Run accession

```bash
$ run=SRR31614351
$ fields=collection_date,country,first_public,location,sample_accession,sample_alias,study_accession
$ curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&query=run_accession=${run}&fields=${fields}"
run_accession   country first_public    location        sample_accession        sample_alias    study_accession collection_date
SRR31614351     Australia:Victoria      2024-12-08              SAMN45172763    AUSMDU00004212  PRJNA856406     2015-11-30
```

---

<!-- TOC --><a name="install-and-use-aspera-to-download-from-enasra"></a>
### Install and use Aspera to download from ENA/SRA

This is all taken from [this amazing tutorial](https://www.biostars.org/p/9528910/).

```
wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0adrj/0/ibm-aspera-connect_4.1.3.93_linux.tar.gz
tar zxvf ibm-aspera-connect_4.1.3.93_linux.tar.gz
bash ibm-aspera-connect_4.1.3.93_linux.sh
$HOME/.aspera/connect/bin/ascp --version
```

Test it out

```
ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR144/004/SRR1448774/SRR1448774.fastq.gz .
```

---

<!-- TOC --><a name="extract-taxon-id-for-list-of-genbank-accessions"></a>
### Extract taxon ID for list of GenBank accessions

(From the master, Wei Shen) extract them from https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

```
zcat nucl_gb.accession2taxid.gz \
    | csvtk grep -t -f accession.version -P accs.txt \
    | csvtk cut -t -f accession.version,taxid
```

or using NCBI tools

```
cat accs.txt
LR789484.1
BC034131.1
BC149740.1
epost -db nucleotide -input accs.txt | esummary | xtract -pattern DocumentSummary -element AccessionVersion,TaxId,Organism
LR789484.1      59560   Phallusia mammillata
BC034131.1      10090   Mus musculus
BC149740.1      9913    Bos taurus
```

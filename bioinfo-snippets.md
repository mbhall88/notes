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

### Convert a VCF file to BED

```
bcftools query -f '%CHROM\t%POS0\t%POS\n' in.vcf
```

---

### Download an assembly, or assemblies, by accession

```
ncbi-genome-download -A GCF_000285655.3 -P -F fasta bacteria
```

`-A` can be a comma separated list of a file to accessions

---

### Change the chromosome name in a VCF

echo -e "old_chr\tnew_chrom" | bcftools annotate --rename-chrs - in.vcf

---

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

---

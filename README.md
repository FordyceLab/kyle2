# Kyle 2.0 - REDI Library Design Utility

<p align="center">
  <img width="200" src="./kyle2.png">
</p>

## How it works

Kyle 2.0 is a design utility for [REDI](http://msb.embopress.org/content/13/2/913.long) libraries. The program works in the following way:

Design library:

1. Define all valid starting point oligos given a set of design parameters.
2. For each starting oligo, find all valid downstream library nodes (where the "window" and primer set are both valid).
3. Repeat until entire region of interested is covered by adjacent, non-overlapping windows.
4. Score each library based on a user-defined scoring function and select best library.

Make mutant oligo:

1. For each desired mutation (either point or insert, scan or user-defined list):
    - Find the codon of interest
    - Change coding sequence to defined codon or add insert **after** codon
    - Write new library oligo to fasta file
    
Mutant oligo names in the fasta file follow the naming convention for point mutations: `{protein name}_{original residue}{position}{mutant residue}_syn:{position of synonymous mutation}_{original codon}->{mutant codon}`. Similarly, oligos with insertions are named: `{protein name}_ins:{position}->{insert}_syn:{position of synonymous mutation}_{original codon}->{mutant codon}`. Examples of formatted versions of these strings include: `pafa_ivtt_Q1A_syn:2_AAA->AAG` and `myc_dhfr_ivtt_ins:145->ACTTGA_syn:145_CAG->CAA` for a point mutation of a PafA IVTT construct and an insertion mutation for a Myc-tagged DHFR IVTT construct.
    
Additionally, all required primers are output to a single CSV file for easy ordering from DNA synthesis companies and follow the following naming convention: `{protein name}_o{oligo ID}_{location (upstream/downstream, encoded as u/d)}_{direction (internal/external, encoded as i/e) }`. This string is appended with `_short` for internal primers that do not include the extended homology tails. An example of the primer name is: `pafa_ivtt_o1_u_i` for an upstream internal primer with hiomology tails amplifies oligo 1 of a library covering the PafA IVTT construct. 

## Installation

Kyle 2.0 requires Python 3 and the cairo backend for SVG production. First install cairo. On a Mac, this can be done with [Homebrew](https://brew.sh/):

```
brew install cairo
```

Then, installation is as easy as:

```
git clone https://github.com/FordyceLab/kyle2.git
cd kyle2
pip3 install .
```

This will make the `kyle2` command line tool publicly available. You can check this with `which kyle2`, which should not return a blank line.

## How to use it

Kyle 2.0 (command line tool `kyle2`) takes a requires set of I/O parameters and has five sets of command line options used to alter its behavior. Of these, the options `--input`, `--output`, `--size`, `--region` and exactly one of \[`--mut_scan`, `--in_scan`, `--mut_file`, or `--insert_file`\] are required:


### I/O flags
- `--input` or `-i` - specifies the input FASTA file containing sequences for oligo design **(required argument)**
- `--output` or `-o` - specifies the output file prefix **(required argument)**
- `--diagram`or `-d` - output an SVG diagram of the best library that was found

### Library type
- `--mut_scan` - mutant residue to scan (single letter amino acid code)
- `--codons` - mutant codon to use in scan (DNA bases)
- `--in_scan` - insert to place at each position (DNA bases)
- `--mut_file` - file specifying mutations to made at specific positions only, 3-column tab-separated no header (position, mutant residue \[1-letter amino acid code\], mutant codon \[DNA bases, optional\])
- `--insert_file` - file specifying insertions to made after specific positions only, 2-column tab-separated no header (position, insert string (DNA bases))

### Library design
- `--size` or `-s` - oligo size **(required argument)**
- `--region` - starting and ending positions of the region of interest (coding sequence) in input fasta file **(required argument)**
- `--primer_length` - primer length range, must provide min and max, default = 16-30
- `--tm` - desired oligo Tm range in degrees C, must provide min and max, default = 54-56
- `--hairpin` - hairpin tolerance in degrees C, default = 10
- `--homodimer` - homodimer tolerance in degrees C, default = 10

### Library scoring options
- `--gc_ends` - GC ends reward for oligo scoring, default = 1
- `--gc_comp` - GC composition penalty for oligo scoring, default = 2
- `--tm_mean` - Tm mean penalty for oligo scoring, default = 1
- `--hairpin_tm` - hairpin Tm suppression reward for oligo scoring, default = 0.1
- `--homodimer_tm` - homodimer Tm suppression reward for oligo scoring, default = 0.1
- `--num_oligos` - set to `min` to only consider libraries with the minimum oligo count, default = `all`

### Tm calculation parameters
- `--dna_conc` - DNA concentration to use for Tm calculation, default = 250
- `--mv_conc` - monovalent cation concentration to use for Tm calculation, default = 50
- `--dv_conc` - divalent cation concentration to use for Tm calculation, default = 0
- `--dntp_conc` - dNTP concentration to use for Tm calculation, default = 0

# Example usage

## Alanine mutational scan

`kyle2 -i pafa_ivtt.fa -o pafa_test_run -d --mut_scan A --size 230 --region 232 1810 --tm 54 55.9 -l 16 30`

## Alanine mutational scan with specified codon

`kyle2 -i pafa_ivtt.fa -o pafa_test_run -d --mut_scan A --codons GCT --size 230 --region 232 1810 --tm 54 55.9 -l 16 30`

## List of specific mutations from a file

`kyle2 -i pafa_ivtt.fa -o pafa_test_run -d --mut_file test_muts.txt --size 230 --region 232 1810 --tm 54 55.9 -l 16 30`

Where the file `test_muts.txt` contains 3 columns corresponding to the position in the construct, the target amino acid residue, and an optional codon to use:

```
1   V
35  M   ATG
225 A   GCA
```

Note that if a codon is specified, it must match the specified amino acid residue.

## Alanine insertion scan

`kyle2 -i pafa_ivtt.fa -o pafa_test_run -d --in_scan GCT --size 230 --region 232 1810 --tm 54 55.9 -l 16 30`

## List of specific insertions from a file

`kyle2 -i pafa_ivtt.fa -o pafa_test_run -d --insert_file test_ins.txt --size 230 --region 232 1810 --tm 54 55.9 -l 16 30`

Where the file `test_ins.txt` contains 2 columns corresponding to the position in the construct **after which** the DNA string will be inserted and the DNA string to insert:

```
24  TTC
67  ATG
122 GCA
```

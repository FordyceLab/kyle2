# Kyle 2.0 - REDI Library Design Utility

<p align="center">
  <img width="200" src="./kyle2.png">
</p>

## How it works

Kyle 2.0 is a design utility for [REDI](http://msb.embopress.org/content/13/2/913.long) libraries. The program works in the following way:

## Installation

Kyle 2.0 requires Python 3.

Installation is as easy as:

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

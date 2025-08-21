
Automated pipeline for protein domain extraction and AlphaFold2 multimer prediction. Converts PDBs to FASTA, extracts variable domains between conserved flanks, and prepares LocalColabFold inputs for structure prediction of domain-partner complexes.

# Protein Domain Extraction & AlphaFold2 Pipeline

![Pipeline](https://img.shields.io/badge/Pipeline-AlphaFold2-blue)
![Python](https://img.shields.io/badge/Python-3.6%2B-green)
![BioPython](https://img.shields.io/badge/BioPython-Optional-orange)
![License](https://img.shields.io/badge/License-MIT-yellow)

Automated pipeline for protein domain extraction and AlphaFold2 multimer prediction. Converts PDBs to FASTA, extracts variable domains between conserved flanking sequences, and prepares LocalColabFold inputs for structure prediction of domain-partner complexes.

## 🎯 What This Pipeline Does

This toolkit streamlines the analysis of designed protein domains by:

1. **Converting protein structures to sequences** with intelligent filename handling
2. **Extracting variable domains** from sequences using conserved flanking regions
3. **Preparing multimer prediction inputs** for LocalColabFold structure prediction

**Perfect for:** Domain grafting projects, protein design validation, binding interface analysis, and structure-function studies of engineered proteins.

## 📋 Requirements

- Python 3.6+
- BioPython (optional, for enhanced PDB parsing)
- LocalColabFold (for structure predictions)

```bash
pip install biopython  # Optional but recommended
```

## 🚀 Quick Start

```bash
# 1. Convert PDB files to FASTA format
python pdb_to_fasta_v2.py

# 2. Extract variable domains using conserved flanks
python domain_extractor.py

# 3. Prepare LocalColabFold multimer inputs
python af2_multimer_prep_type2_local.py

# 4. Run structure predictions with LocalColabFold
--> https://github.com/YoshitakaMo/localcolabfold
```

## 📖 Detailed Script Documentation

### 1. `pdb_to_fasta_v2.py`
**Purpose:** Converts PDB files to FASTA format with intelligent filename processing

**Configuration (edit top of script):**
```python
PDB_DIRECTORY = "/path/to/your/pdb/files"  # Directory containing PDB files
CHAIN_ID = "A"                             # Chain to extract (default: A)
```

**Key Features:**
- Handles complex filenames: `6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.pdb`
- Creates simplified headers: `>6_dir6_n1-3_20250705_33_7`
- Outputs individual `.fasta` files + combined `all_sequences.txt`
- Works with or without BioPython

**Input:** Directory of PDB files  
**Output:** Individual FASTA files + `all_sequences.txt`

**Example Usage:**
```bash
# Basic usage (edit PDB_DIRECTORY in script first)
python pdb_to_fasta_v2.py

# Output files created:
#   individual_file.fasta (for each PDB)
#   all_sequences.txt (combined file)
```

---

### 2. `domain_extractor.py`
**Purpose:** Extracts variable domain regions between conserved flanking sequences

**Configuration (edit top of script):**
```python
INPUT_FASTA = "all_sequences.txt"           # Input from previous script
OUTPUT_FASTA = "extracted_domains.txt"      # Output file name
N_TERM_FLANK = "VYTEDEWQKEWNELIKLASSEP"     # N-terminal conserved sequence
C_TERM_FLANK = "EPVYESLEEFHVFVLAHVLRRP"     # C-terminal conserved sequence  
KEEP_RESIDUES = 5                           # Residues to keep from each flank
```

**Key Features:**
- Finds variable regions between conserved sequences
- Keeps specified residues from constant regions for structural context
- Validates flanking sequences are present
- Detailed extraction preview and logging

**Extraction Logic:**
```
Full sequence: ----[N_FLANK]===VARIABLE_DOMAIN===[C_FLANK]----
Extracted:          [last 5]===VARIABLE_DOMAIN===[first 5]
```

**Input:** FASTA file with full sequences  
**Output:** FASTA file with extracted domains

**Example Usage:**
```bash
# Basic usage (configure flanking sequences first)
python domain_extractor.py

# The script will show preview:
# N-terminal: ...VYTEDEWQKEWNELIKLASSEP
#             Keep last 5: ASSEP
# C-terminal: EPVYESLEEFHVFVLAHVLRRP...
#             Keep first 5: EPVYE
```

---

### 3. `af2_multimer_prep_type2_local.py`
**Purpose:** Prepares LocalColabFold multimer prediction inputs by combining domains with partner proteins

**Configuration (edit top of script):**
```python
INPUT_DOMAINS = "extracted_domains.txt"           # Domains from previous script
PARTNER_SEQUENCE = """MQIFVKTLTGKT..."""         # Full partner protein sequence
OUTPUT_DIRECTORY = "colabfold_multimer_inputs"    # Directory for individual files
COMBINED_OUTPUT = "all_colabfold_multimer_inputs.fasta"  # Combined file
SEQUENCE_ORDER = "partner_first"                  # "partner_first" or "domain_first"
```

**Key Features:**
- Creates colon-separated multimer format for LocalColabFold
- Combines each domain with partner protein
- Supports both individual and batch file creation
- Preserves domain names for tracking variants

**Output Format:**
```
>domain_name
PARTNER_PROTEIN_SEQUENCE:DOMAIN_SEQUENCE
```

**Input:** Extracted domains + partner sequence  
**Output:** LocalColabFold-ready multimer inputs

**Example Usage:**
```bash
# Configure partner sequence in script, then run:
python af2_multimer_prep_type2_local.py

# Creates files ready for:
colabfold_batch all_colabfold_multimer_inputs.fasta predictions/
```

## 📁 Complete Workflow Example

### Input Files:
```
my_designs/
├── design1.pdb
├── design2.pdb
├── design3.pdb
└── ...
```

### Step 1: PDB to FASTA Conversion
```bash
# Edit pdb_to_fasta_v2.py:
PDB_DIRECTORY = "my_designs/"

python pdb_to_fasta_v2.py
```

**Creates:**
```
my_designs/
├── design1.fasta
├── design2.fasta  
├── design3.fasta
├── all_sequences.txt  ← Input for next step
└── ...
```

### Step 2: Domain Extraction
```bash
# Edit domain_extractor.py:
INPUT_FASTA = "my_designs/all_sequences.txt"
N_TERM_FLANK = "VYTEDEWQKEWNELIKLASSEP"  # Your conserved N-terminal
C_TERM_FLANK = "EPVYESLEEFHVFVLAHVLRRP"  # Your conserved C-terminal

python domain_extractor.py
```

**Creates:**
```
extracted_domains.txt  ← Input for next step
```

### Step 3: Multimer Input Preparation
```bash
# Edit af2_multimer_prep_type2_local.py:
PARTNER_SEQUENCE = """MQIFVKTLTGKTITLE..."""  # Your partner protein

python af2_multimer_prep_type2_local.py
```

**Creates:**
```
colabfold_multimer_inputs/
├── domain1.fasta
├── domain2.fasta
├── domain3.fasta
└── ...
all_colabfold_multimer_inputs.fasta  ← Ready for prediction
```

### Step 4: Structure Prediction
```bash
# Run LocalColabFold predictions
--> https://github.com/YoshitakaMo/localcolabfold
```

## ⚙️ Configuration Guide

### Domain Extractor Settings

**Finding Conserved Flanks:**
1. Align your full-length sequences
2. Identify conserved regions before/after variable domain
3. Choose 15-25 residue flanking sequences
4. Set `KEEP_RESIDUES` to maintain structural context (usually 3-7)

**Example Configuration:**
```python
# For immunoglobulin-like domains:
N_TERM_FLANK = "VYTEDEWQKEWNELIKLASSEP"  # 21 residues
C_TERM_FLANK = "EPVYESLEEFHVFVLAHVLRRP"  # 21 residues
KEEP_RESIDUES = 5  # Keep 5 from each end
```

### Multimer Prep Settings

**Sequence Ordering:**
- `"partner_first"`: `PARTNER:DOMAIN` (recommended for most cases)
- `"domain_first"`: `DOMAIN:PARTNER` (if domain is larger/more structured)

**File Organization:**
- Individual files: Better for parallel processing
- Combined file: Easier for single batch runs

## 🔧 File Flow Diagram

```
PDB Files ──→ [Script 1] ──→ all_sequences.txt
                                      ↓
Conserved Flanks ──→ [Script 2] ──→ extracted_domains.txt
                                            ↓
Partner Sequence ──→ [Script 3] ──→ multimer_inputs.fasta
                                            ↓
                              LocalColabFold ──→ Predicted Structures
```

## 🔬 Scientific Applications

**Perfect for:**
- **Domain grafting validation** - Check if designed domains maintain binding
- **Binding interface analysis** - Study domain-partner interactions
- **Design optimization** - Compare multiple domain variants
- **Conformational studies** - Understand domain flexibility in context

**Typical Use Cases:**
- Antibody domain engineering
- Enzyme active site grafting  
- Protein-protein interaction design
- Modular protein assembly validation

## 🐛 Troubleshooting

### Common Issues

**Script 1 - No sequences extracted:**
- Check `PDB_DIRECTORY` path
- Verify `CHAIN_ID` exists in PDB files
- Check file permissions

**Script 2 - No domains found:**
- Verify flanking sequences are present in input
- Check sequence case sensitivity
- Try shorter flanking sequences
- Increase `KEEP_RESIDUES` if domains seem too short

**Script 3 - Format errors:**
- Ensure partner sequence is valid amino acids
- Check for special characters in domain names
- Verify input domain file exists

**LocalColabFold issues:**
- Use `--model-type auto` for automatic multimer detection
- Add `--amber` for structure relaxation
- Check sequence lengths (very long sequences may fail)

### Debugging Tips

**Validate extractions:**
```bash
# Check extraction results
grep ">" extracted_domains.txt | head -5
grep -v ">" extracted_domains.txt | head -1 | wc -c
```

**Test multimer format:**
```bash
# Verify colon-separated format
head -2 all_colabfold_multimer_inputs.fasta
```

## 📊 Output Quality Checks

**After Script 1:**
- Verify all PDB files were processed
- Check sequence lengths are reasonable
- Confirm headers are properly formatted

**After Script 2:**  
- Ensure extracted domains start/end with expected residues
- Check domain lengths are consistent with expectations
- Verify no domains were lost due to missing flanks

**After Script 3:**
- Confirm colon-separated format: `SEQUENCE1:SEQUENCE2`
- Check total sequence lengths per multimer
- Verify domain names are preserved in headers

## 📧 Contact

e.leen@leeds.ac.uk

## 📄 License

MIT License - see LICENSE file for details

## 🔗 Related Tools

- [LocalColabFold](https://github.com/YoshitakaMo/localcolabfold) - Local AlphaFold2 predictions
- [ColabFold](https://colab.research.google.com/github/deepmind/ColabFold/) - Google Colab notebooks
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) - Structure visualization
- [PyMOL](https://pymol.org/) - Molecular graphics

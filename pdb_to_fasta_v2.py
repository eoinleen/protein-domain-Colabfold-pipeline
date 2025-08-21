# PDB to FASTA Converter for Red Hat Linux - Modified Version
# =================================================
# This script extracts amino acid sequences from PDB files
# and saves them in FASTA format. Adapted for complex filenames.
#
# Features:
# - Handles complex filename patterns like: 6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.pdb
# - Creates simplified FASTA headers: >6_dir6_n1-3_20250705_33_7
# - Outputs .fasta files with same base name as input
# - Creates a combined 'all_sequences.txt' FASTA file
#
# Requirements:
# - Python 3.6+ with Biopython installed (optional)
# - Run: pip install biopython (or use conda/yum as needed)
#
# Usage:
# - Set 'PDB_DIRECTORY' to the path containing PDB files
# - Set 'CHAIN_ID' to the chain you want to extract (default: A)
# - Run: python pdb_to_fasta.py
# =================================================

# ===== SETTINGS - EDIT THESE =====
PDB_DIRECTORY = "../all_extracted_pdb"  # Update this path
CHAIN_ID = "A"
# ================================

import os
import re
import glob
import sys

# Check if biopython is available (optional - script works without it)
try:
    from Bio import PDB
    BIOPYTHON_AVAILABLE = True
    print("Biopython detected - enhanced PDB parsing available")
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Biopython not found - using basic PDB parsing")
    print("To install: pip install biopython")

# Dictionary to convert three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence_from_pdb(pdb_file, chain_id):
    """Extract the amino acid sequence from a specified chain in a PDB file."""
    sequence = ""
    seen_residues = set()

    try:
        with open(pdb_file, 'r', encoding="utf-8") as f:
            for line in f:
                if line.startswith('ATOM') and len(line) > 21 and line[21] == chain_id:
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()

                    if res_num not in seen_residues:
                        if res_name in three_to_one:
                            sequence += three_to_one[res_name]
                        seen_residues.add(res_num)
    except Exception as e:
        print(f"Error reading {pdb_file}: {e}")
        return ""

    return sequence


def create_fasta_header(filename):
    """Create a simplified FASTA header from complex filename.
    
    Example: 6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.pdb
    Returns: >6_dir6_n1-3_20250705_33_7
    """
    # Remove .pdb extension
    base_name = os.path.splitext(filename)[0]
    
    # Split by underscores
    parts = base_name.split('_')
    
    header_parts = []
    
    for part in parts:
        # Skip common suffixes
        if part.lower() in ['af2pred', 'dldesign']:
            continue
            
        # Simplify 'noise1-3' to 'n1-3'
        if part.startswith('noise'):
            simplified = 'n' + part[5:]  # Remove 'noise', keep the rest
            header_parts.append(simplified)
        else:
            header_parts.append(part)
    
    # Join with underscores and add > prefix
    header = '>' + '_'.join(header_parts)
    return header


def create_output_filename(pdb_filename):
    """Create output .fasta filename from input .pdb filename.
    
    Example: 6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.pdb
    Returns: 6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.fasta
    """
    base_name = os.path.splitext(pdb_filename)[0]
    return base_name + ".fasta"


def convert_pdb_to_fasta(pdb_dir, chain_id):
    """Convert all PDB files in the directory to FASTA format."""
    if not os.path.isdir(pdb_dir):
        print(f"Error: Directory not found: {pdb_dir}")
        print("Please update the PDB_DIRECTORY variable with the correct path.")
        return False

    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}")
        return False

    print(f"Processing {len(pdb_files)} PDB files...")
    all_fasta_entries = []
    processed_count = 0

    for pdb_file in pdb_files:
        filename = os.path.basename(pdb_file)
        
        # Create FASTA header and output filename
        fasta_header = create_fasta_header(filename)
        output_filename = create_output_filename(filename)
        
        # Extract sequence
        sequence = extract_sequence_from_pdb(pdb_file, chain_id)

        if not sequence:
            print(f"Warning: No sequence found for {filename}")
            continue

        # Format FASTA entry with 60 characters per line
        fasta_entry = fasta_header + "\n" + "\n".join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
        all_fasta_entries.append(fasta_entry)

        # Write individual FASTA file
        output_file = os.path.join(pdb_dir, output_filename)
        try:
            with open(output_file, 'w', encoding="utf-8") as f:
                f.write(fasta_entry)
            print(f"Processed {filename} -> {output_filename}")
            print(f"  Header: {fasta_header}")
            processed_count += 1
        except Exception as e:
            print(f"Error writing {output_file}: {e}")

    # Create combined file
    if processed_count > 0:
        all_sequences_file = os.path.join(pdb_dir, "all_sequences.txt")
        try:
            with open(all_sequences_file, 'w', encoding="utf-8") as f:
                f.write("\n".join(all_fasta_entries))
            print(f"\nCreated combined FASTA file: {all_sequences_file}")
        except Exception as e:
            print(f"Error creating combined file: {e}")

    print(f"\nProcessed {processed_count} PDB files.")
    return True


def main():
    """Main execution function."""
    print("PDB to FASTA Converter for Red Hat Linux - Modified Version")
    print("="*60)

    if PDB_DIRECTORY == "/path/to/your/pdb/files":
        print("WARNING: Please update the PDB_DIRECTORY variable with your actual path!")
        print("Edit the script and change the PDB_DIRECTORY setting at the top.")
        return

    print(f"PDB Directory: {PDB_DIRECTORY}")
    print(f"Chain ID: {CHAIN_ID}")
    print()

    # Run the conversion
    success = convert_pdb_to_fasta(PDB_DIRECTORY, CHAIN_ID)

    if success:
        print("\nConversion completed successfully!")
        print("Example transformations:")
        print("  Input:  6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.pdb")
        print("  Output: 6_dir6_noise1-3_20250705_33_dldesign_7_af2pred.fasta")
        print("  Header: >6_dir6_n1-3_20250705_33_7")
    else:
        print("\nConversion failed. Please check the directory path and file permissions.")


if __name__ == "__main__":
    main()

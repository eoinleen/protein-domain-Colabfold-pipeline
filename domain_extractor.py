# Variable Domain Extractor
# =================================================
# This script extracts variable domain regions from protein sequences
# based on conserved flanking sequences. It reads FASTA files and
# extracts the region between two specified sequences plus a defined
# number of residues from the constant regions.
#
# Usage:
# - Set input and output file paths
# - Define N-terminal and C-terminal flanking sequences
# - Set number of residues to keep from constant regions
# - Run: python domain_extractor.py
# =================================================

# ===== SETTINGS - EDIT THESE =====
INPUT_FASTA = "all_sequences.txt"  # Input FASTA file from previous script
OUTPUT_FASTA = "extracted_domains.txt"  # Output file for extracted sequences

# Flanking sequences (without the overlap residues)
N_TERM_FLANK = "VYTEDEWQKEWNELIKLASSEP"  # N-terminal constant region
C_TERM_FLANK = "EPVYESLEEFHVFVLAHVLRRP"  # C-terminal constant region

# Number of residues to keep from each constant region
KEEP_RESIDUES = 5
# ================================

import os
import sys

def read_fasta(filename):
    """Read FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_sequence = ""
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence if it exists
                    if current_header is not None:
                        sequences.append((current_header, current_sequence))
                    # Start new sequence
                    current_header = line
                    current_sequence = ""
                else:
                    current_sequence += line
            
            # Save the last sequence
            if current_header is not None:
                sequences.append((current_header, current_sequence))
                
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return []
    except Exception as e:
        print(f"Error reading file '{filename}': {e}")
        return []
    
    return sequences


def extract_domain(sequence, n_term_flank, c_term_flank, keep_residues):
    """Extract variable domain between flanking sequences.
    
    Args:
        sequence: Full protein sequence
        n_term_flank: N-terminal flanking sequence
        c_term_flank: C-terminal flanking sequence  
        keep_residues: Number of residues to keep from each flank
        
    Returns:
        Extracted domain sequence or None if flanks not found
    """
    # Find positions of flanking sequences
    n_term_pos = sequence.find(n_term_flank)
    c_term_pos = sequence.find(c_term_flank)
    
    if n_term_pos == -1:
        print(f"Warning: N-terminal flank '{n_term_flank}' not found")
        return None
        
    if c_term_pos == -1:
        print(f"Warning: C-terminal flank '{c_term_flank}' not found")
        return None
    
    if c_term_pos <= n_term_pos:
        print(f"Warning: C-terminal flank appears before N-terminal flank")
        return None
    
    # Calculate extraction positions
    # Start: end of N-terminal flank minus keep_residues
    start_pos = n_term_pos + len(n_term_flank) - keep_residues
    
    # End: start of C-terminal flank plus keep_residues  
    end_pos = c_term_pos + keep_residues
    
    # Extract the domain
    extracted = sequence[start_pos:end_pos]
    
    return extracted


def format_fasta_sequence(sequence, line_length=60):
    """Format sequence with specified line length."""
    return '\n'.join([sequence[i:i+line_length] for i in range(0, len(sequence), line_length)])


def extract_domains_from_fasta(input_file, output_file, n_term_flank, c_term_flank, keep_residues):
    """Process FASTA file and extract domains from all sequences."""
    
    # Read input sequences
    sequences = read_fasta(input_file)
    if not sequences:
        return False
    
    print(f"Processing {len(sequences)} sequences...")
    print(f"N-terminal flank: {n_term_flank}")
    print(f"C-terminal flank: {c_term_flank}")
    print(f"Keeping {keep_residues} residues from each flank")
    print()
    
    extracted_sequences = []
    successful_extractions = 0
    
    for header, sequence in sequences:
        print(f"Processing {header}...")
        
        # Extract domain
        extracted = extract_domain(sequence, n_term_flank, c_term_flank, keep_residues)
        
        if extracted is not None:
            extracted_sequences.append((header, extracted))
            successful_extractions += 1
            print(f"  Extracted {len(extracted)} residues: {extracted[:20]}{'...' if len(extracted) > 20 else ''}")
        else:
            print(f"  Failed to extract domain")
        print()
    
    # Write output file
    if extracted_sequences:
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                for header, sequence in extracted_sequences:
                    f.write(header + '\n')
                    f.write(format_fasta_sequence(sequence) + '\n')
            
            print(f"Successfully extracted domains from {successful_extractions}/{len(sequences)} sequences")
            print(f"Output written to: {output_file}")
            return True
            
        except Exception as e:
            print(f"Error writing output file: {e}")
            return False
    else:
        print("No domains were successfully extracted.")
        return False


def validate_settings():
    """Validate the settings before processing."""
    issues = []
    
    if not os.path.exists(INPUT_FASTA):
        issues.append(f"Input file '{INPUT_FASTA}' does not exist")
    
    if len(N_TERM_FLANK) < KEEP_RESIDUES:
        issues.append(f"N-terminal flank ({len(N_TERM_FLANK)} residues) is shorter than KEEP_RESIDUES ({KEEP_RESIDUES})")
    
    if len(C_TERM_FLANK) < KEEP_RESIDUES:
        issues.append(f"C-terminal flank ({len(C_TERM_FLANK)} residues) is shorter than KEEP_RESIDUES ({KEEP_RESIDUES})")
    
    if KEEP_RESIDUES < 0:
        issues.append("KEEP_RESIDUES must be non-negative")
    
    return issues


def main():
    """Main execution function."""
    print("Variable Domain Extractor")
    print("=" * 50)
    print()
    
    # Validate settings
    issues = validate_settings()
    if issues:
        print("Configuration issues found:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nPlease fix these issues and try again.")
        return
    
    # Show extraction preview
    print("Extraction Preview:")
    print(f"N-terminal: ...{N_TERM_FLANK}")
    print(f"            Keep last {KEEP_RESIDUES}: {N_TERM_FLANK[-KEEP_RESIDUES:]}")
    print(f"C-terminal: {C_TERM_FLANK}...")
    print(f"            Keep first {KEEP_RESIDUES}: {C_TERM_FLANK[:KEEP_RESIDUES]}")
    print(f"Extracted sequence will start with: {N_TERM_FLANK[-KEEP_RESIDUES:]}")
    print(f"Extracted sequence will end with: {C_TERM_FLANK[:KEEP_RESIDUES]}")
    print()
    
    # Process files
    success = extract_domains_from_fasta(
        INPUT_FASTA, 
        OUTPUT_FASTA, 
        N_TERM_FLANK, 
        C_TERM_FLANK, 
        KEEP_RESIDUES
    )
    
    if success:
        print("\nDomain extraction completed successfully!")
    else:
        print("\nDomain extraction failed.")


if __name__ == "__main__":
    main()

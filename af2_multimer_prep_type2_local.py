# AF2 Multimer Input Preparation Script - Type 2 (LocalColabFold Format)
# =================================================================================
# DETAILED OVERVIEW:
# This script prepares input files for LocalColabFold multimer predictions by combining
# extracted variable domain sequences with a full-length partner protein using the
# colon-separated format required by LocalColabFold/ColabFold systems.
#
# WORKFLOW PROCESS:
# 1. INPUT: Reads extracted domain sequences from domain extraction script
#    - Example input domains: >1_dir1_n0-8_20250705_1_26
#                            ASSEPRGEERFERETIKSLERILESLK...EPVYE
#
# 2. PARTNER: Takes user-provided full-length protein sequence
#    - Example: Full scaffolded protein (constant + variable regions)
#
# 3. COMBINATION: Creates multimer complexes using colon-separated format
#    - Domain 1 + Full Protein = Complex 1 (FULL_PROTEIN:DOMAIN_1)
#    - Domain 2 + Full Protein = Complex 2 (FULL_PROTEIN:DOMAIN_2)
#    - Domain N + Full Protein = Complex N (FULL_PROTEIN:DOMAIN_N)
#
# 4. OUTPUT FORMAT: LocalColabFold colon-separated format (single header, colon between chains)
#    - Each complex saved as individual .fasta file
#    - All complexes optionally combined in single batch file
#    - Domain name preserved as header (not "full_protein")
#
# DETAILED INPUT/OUTPUT EXAMPLE:
# 
# INPUT (from extracted_domains.txt):
# >1_dir1_n0-8_20250705_1_26
# ASSEPRGEERFERETIKSLERILESLKRLKKIAEKEGKKEKAEKYEKEAAEKEKELAAKK
# AEFAKVAPLDTEPVYE
# >2_dir2_n1-5_20250705_2_15
# ASSEPRGAGEPEALEEAMAALAAEAGLDPAEAAARMRELYAAGRDVEALREWLVATFGDE
# AAAEELAALVRRALEEPVYE
#
# OUTPUT - Individual Files:
# File: 1_dir1_n0-8_20250705_1_26.fasta
# >1_dir1_n0-8_20250705_1_26
# MEMPICAFQLPDLTVYNEDFRSFIERDLIEQSMLVALEQAGRLNWWVSVDPTSQRLLPLA
# TTGDGNCLLHAASLGMWGFHDRDLMLRKALYALMEKGVEKEALKRRWRWQQTQQNKESGL
# VYTEDEWQKEWNELIKLASSEPRGPEVLEKAKEELEKKIEELKAEDSEEAEKEAEAAKLM
# LEALETKEEPVYESLEEFHVFVLAHVLRRPIVVVADTMLRGPGGEAFAPIPFGGIYLPLE
# VPASQCHRSPLVLAYDQAHFSALVSMEQKENTKEQAVIPLTDSEYKLLPLHFAVDPGKGW
# EWGKDDSDNVRLASVILSLEVKLHLLHSYMNVKWIP:ASSEPRGEERFERETIKSLERIL
# ESLKRLKKIAEKEGKKEKAEKYEKEAAEKEKELAAKKAEFAKVAPLDTEPVYE
#
# File: 2_dir2_n1-5_20250705_2_15.fasta
# >2_dir2_n1-5_20250705_2_15
# MEMPICAFQLPDLTVYNEDFRSFIERDLIEQSMLVALEQAGRLNWWVSVDPTSQRLLPLA
# TTGDGNCLLHAASLGMWGFHDRDLMLRKALYALMEKGVEKEALKRRWRWQQTQQNKESGL
# VYTEDEWQKEWNELIKLASSEPRGPEVLEKAKEELEKKIEELKAEDSEEAEKEAEAAKLM
# LEALETKEEPVYESLEEFHVFVLAHVLRRPIVVVADTMLRGPGGEAFAPIPFGGIYLPLE
# VPASQCHRSPLVLAYDQAHFSALVSMEQKENTKEQAVIPLTDSEYKLLPLHFAVDPGKGW
# EWGKDDSDNVRLASVILSLEVKLHLLHSYMNVKWIP:ASSEPRGAGEPEALEEAMAALA
# AEGLDPAEAAARMRELYAAGRDVEALREWLVATFGDEAAAEELAALVRRALEEPVYE
#
# OUTPUT - Combined File (all_multimer_inputs.fasta):
# >1_dir1_n0-8_20250705_1_26
# MEMPICAFQL...VKWIP:ASSEPRGEERFERE...EPVYE
# >2_dir2_n1-5_20250705_2_15
# MEMPICAFQL...VKWIP:ASSEPRGAGEPEAL...EPVYE
# [... continues for all domains]
#
# SEQUENCE TRANSFORMATIONS:
# - Each extracted domain becomes second chain in colon-separated complex
# - Partner protein sequence is repeated as first chain for every complex
# - Original domain headers are preserved for tracking design variants
# - Single line format (no 60-character wrapping) for colon-separated sequences
# - Colon (:) separates the two protein chains within single sequence line
#
# SCIENTIFIC PURPOSE:
# This creates input files to predict how each designed domain variant
# interacts with the full scaffolded protein using LocalColabFold's
# efficient multimer prediction pipeline.
#
# COMPATIBLE SYSTEMS:
# - LocalColabFold (YoshitakaMo/localcolabfold)
# - ColabFold (Google Colab notebooks)
# - Any system expecting colon-separated multimer format
#
# NOT COMPATIBLE WITH:
# - Native AlphaFold2 installations (requires separate headers - see Type 1 script)
# - Docker-based AlphaFold2 (requires separate headers - see Type 1 script)
#
# USAGE:
# 1. Set INPUT_DOMAINS to output file from domain extraction script
# 2. Paste full protein sequence into PARTNER_SEQUENCE variable
# 3. Configure output settings (individual files, combined file, naming)
# 4. Run: python af2_multimer_prep_type2.py
# 5. Use output files with LocalColabFold: colabfold_batch input_dir/ output_dir/
#
# BATCH PROCESSING EXAMPLE:
# After running this script, use with multirun.sh style batch processing:
# for fastaname in domain1 domain2 domain3 ; do
#     colabfold_batch $fastaname.fasta outputs_$fastaname/
# done
# =================================================================================

# ===== SETTINGS - EDIT THESE =====
# Input file from domain extractor
INPUT_DOMAINS = "extracted_domains.txt"

# Partner sequence to combine with domains
PARTNER_SEQUENCE = """MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"""

# Output settings
OUTPUT_DIRECTORY = "colabfold_multimer_inputs"  # Directory for individual files
COMBINED_OUTPUT = "all_colabfold_multimer_inputs.fasta"  # Combined file
CREATE_INDIVIDUAL_FILES = True  # Create separate files for each pair
CREATE_COMBINED_FILE = True     # Create one file with all pairs

# Sequence ordering: "domain_first" or "partner_first"
SEQUENCE_ORDER = "partner_first"

# File naming format for individual files
# Options: "domain_name", "numbered", "custom"
NAMING_FORMAT = "domain_name"  # Uses domain header for filename
# ================================

import os
import sys
import re

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
                    current_header = line[1:]  # Remove '>' prefix
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


def clean_sequence(sequence):
    """Remove any whitespace or formatting from sequence."""
    return ''.join(sequence.split())


def create_safe_filename(name, extension=".fasta"):
    """Create a safe filename from a sequence name."""
    # Remove problematic characters and replace with underscores
    safe_name = re.sub(r'[<>:"/\\|?*]', '_', name)
    safe_name = re.sub(r'[^\w\-_.]', '_', safe_name)
    # Remove multiple consecutive underscores
    safe_name = re.sub(r'_{2,}', '_', safe_name)
    # Remove leading/trailing underscores
    safe_name = safe_name.strip('_')
    
    return safe_name + extension


def create_colabfold_multimer_input(domain_header, domain_sequence, partner_sequence, sequence_order):
    """Create a LocalColabFold multimer input with colon-separated sequences."""
    
    if sequence_order == "domain_first":
        colabfold_sequence = f"{domain_sequence}:{partner_sequence}"
    else:  # partner_first
        colabfold_sequence = f"{partner_sequence}:{domain_sequence}"
    
    fasta_content = f">{domain_header}\n{colabfold_sequence}\n"
    
    return fasta_content


def prepare_colabfold_multimer_inputs(input_file, partner_seq, output_dir, combined_file, 
                                    create_individual, create_combined, sequence_order, naming_format):
    """Process domain sequences and create LocalColabFold multimer inputs."""
    
    # Read domain sequences
    domains = read_fasta(input_file)
    if not domains:
        print("No domain sequences found. Exiting.")
        return False
    
    # Clean partner sequence
    clean_partner_seq = clean_sequence(partner_seq)
    
    print(f"Processing {len(domains)} domain sequences...")
    print(f"Partner sequence: {len(clean_partner_seq)} residues")
    print(f"Sequence order: {sequence_order}")
    print(f"Output format: LocalColabFold colon-separated")
    print()
    
    # Create output directory if needed
    if create_individual:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory: {output_dir}")
    
    all_multimer_inputs = []
    successful_count = 0
    
    for i, (domain_header, domain_sequence) in enumerate(domains, 1):
        clean_domain_seq = clean_sequence(domain_sequence)
        
        print(f"Processing {i}/{len(domains)}: {domain_header}")
        print(f"  Domain length: {len(clean_domain_seq)} residues")
        
        # Create LocalColabFold multimer input
        colabfold_fasta = create_colabfold_multimer_input(
            domain_header, clean_domain_seq, 
            clean_partner_seq, sequence_order
        )
        
        # Add to combined list
        all_multimer_inputs.append(colabfold_fasta.strip())
        
        # Create individual file if requested
        if create_individual:
            if naming_format == "domain_name":
                filename = create_safe_filename(domain_header)
            elif naming_format == "numbered":
                filename = f"colabfold_multimer_{i:03d}.fasta"
            else:  # custom - you can modify this
                filename = f"{domain_header}_colabfold.fasta"
                filename = create_safe_filename(filename)
            
            filepath = os.path.join(output_dir, filename)
            
            try:
                with open(filepath, 'w', encoding='utf-8') as f:
                    f.write(colabfold_fasta)
                print(f"  Created: {filename}")
                successful_count += 1
            except Exception as e:
                print(f"  Error creating {filename}: {e}")
        else:
            successful_count += 1
    
    # Create combined file if requested
    if create_combined and all_multimer_inputs:
        try:
            with open(combined_file, 'w', encoding='utf-8') as f:
                f.write('\n'.join(all_multimer_inputs) + '\n')
            print(f"\nCreated combined file: {combined_file}")
            print(f"Contains {len(all_multimer_inputs)} multimer inputs")
        except Exception as e:
            print(f"Error creating combined file: {e}")
            return False
    
    print(f"\nSuccessfully processed {successful_count}/{len(domains)} domain sequences")
    return True


def validate_settings():
    """Validate settings before processing."""
    issues = []
    
    if not os.path.exists(INPUT_DOMAINS):
        issues.append(f"Input file '{INPUT_DOMAINS}' does not exist")
    
    if not PARTNER_SEQUENCE.strip():
        issues.append("PARTNER_SEQUENCE is empty")
    
    if SEQUENCE_ORDER not in ["domain_first", "partner_first"]:
        issues.append("SEQUENCE_ORDER must be 'domain_first' or 'partner_first'")
    
    if NAMING_FORMAT not in ["domain_name", "numbered", "custom"]:
        issues.append("NAMING_FORMAT must be 'domain_name', 'numbered', or 'custom'")
    
    if not CREATE_INDIVIDUAL_FILES and not CREATE_COMBINED_FILE:
        issues.append("At least one of CREATE_INDIVIDUAL_FILES or CREATE_COMBINED_FILE must be True")
    
    return issues


def main():
    """Main execution function."""
    print("AF2 Multimer Input Preparation Script - Type 2 (LocalColabFold Format)")
    print("=" * 70)
    print()
    
    # Validate settings
    issues = validate_settings()
    if issues:
        print("Configuration issues found:")
        for issue in issues:
            print(f"  - {issue}")
        print("\nPlease fix these issues and try again.")
        return
    
    # Show configuration
    print("Configuration:")
    print(f"  Input domains: {INPUT_DOMAINS}")
    print(f"  Partner length: {len(clean_sequence(PARTNER_SEQUENCE))} residues")
    print(f"  Sequence order: {SEQUENCE_ORDER}")
    print(f"  Create individual files: {CREATE_INDIVIDUAL_FILES}")
    if CREATE_INDIVIDUAL_FILES:
        print(f"    Output directory: {OUTPUT_DIRECTORY}")
        print(f"    Naming format: {NAMING_FORMAT}")
    print(f"  Create combined file: {CREATE_COMBINED_FILE}")
    if CREATE_COMBINED_FILE:
        print(f"    Combined file: {COMBINED_OUTPUT}")
    print()
    
    # Process files
    success = prepare_colabfold_multimer_inputs(
        INPUT_DOMAINS, PARTNER_SEQUENCE,
        OUTPUT_DIRECTORY, COMBINED_OUTPUT,
        CREATE_INDIVIDUAL_FILES, CREATE_COMBINED_FILE,
        SEQUENCE_ORDER, NAMING_FORMAT
    )
    
    if success:
        print("\nLocalColabFold multimer input preparation completed successfully!")
        print("\nNext steps for LocalColabFold:")
        if CREATE_INDIVIDUAL_FILES:
            print(f"  - Process individual files from '{OUTPUT_DIRECTORY}/' directory")
            print("  - Example batch command:")
            print("    for fastaname in domain1 domain2 domain3 ; do")
            print("        colabfold_batch $fastaname.fasta outputs_$fastaname/")
            print("    done")
        if CREATE_COMBINED_FILE:
            print(f"  - Or process combined file: colabfold_batch {COMBINED_OUTPUT} output_dir/")
        print("  - Use --model-type auto (default) for automatic multimer detection")
        print("  - Add --amber for structure relaxation (recommended)")
    else:
        print("\nLocalColabFold multimer input preparation failed.")


if __name__ == "__main__":
    main()

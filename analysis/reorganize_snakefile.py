#!/usr/bin/env python3
"""
Script to reorganize Snakefile rules in logical data flow order.
"""

import re
import sys

# Define the desired order of rules
RULE_ORDER = [
    # Keep preamble (imports, config, etc.)
    "PREAMBLE",
    
    # Main rule
    "all",
    
    # STAGE 1: Input processing and counting
    "count_reads",
    "calculate_pool_statistics",
    
    # STAGE 2: Undetermined recovery and barcode correction
    "check_undetermined_barcodes",
    "create_undetermined_fastq", 
    "barcode_recovery",
    "merge_all_fastqs",
    
    # STAGE 3: Reference preparation
    "kite_index",
    "generate_combined_whitelist",
    
    # STAGE 4: Alignment
    "kallisto_gex",
    "kallisto_guide",
    "kallisto_gex_subsampled",
    "kallisto_guide_subsampled",
    
    # STAGE 5: Post-alignment processing
    "inspect_bus_files",
    "filter_and_annotate_sublibrary",
    "calculate_read_statistics",
    
    # STAGE 6: QC and analysis
    "fastp_qc",
    "cell_calling_analysis",
    "cell_calling_plots",
    "calculate_qc_metrics_stratified",
    "umi_saturation_analysis",
    "umi_saturation_analysis_guide",
    
    # STAGE 7: Consolidation and visualization
    "consolidate_qc_metrics",
    "visualize_consolidated_qc",
    "process_pool_metrics",
    
    # STAGE 8: Final outputs
    "generate_qc_report",
    "combine_sublibraries",
]

# Section headers to insert
SECTION_HEADERS = {
    "count_reads": "\n# =============================================================================\n# STAGE 1: INPUT PROCESSING AND COUNTING\n# =============================================================================\n",
    "check_undetermined_barcodes": "\n# =============================================================================\n# STAGE 2: UNDETERMINED RECOVERY AND BARCODE CORRECTION\n# =============================================================================\n",
    "kite_index": "\n# =============================================================================\n# STAGE 3: REFERENCE PREPARATION\n# =============================================================================\n",
    "kallisto_gex": "\n# =============================================================================\n# STAGE 4: ALIGNMENT (KALLISTO/BUSTOOLS)\n# =============================================================================\n",
    "inspect_bus_files": "\n# =============================================================================\n# STAGE 5: POST-ALIGNMENT PROCESSING\n# =============================================================================\n",
    "fastp_qc": "\n# =============================================================================\n# STAGE 6: QC AND ANALYSIS\n# =============================================================================\n",
    "consolidate_qc_metrics": "\n# =============================================================================\n# STAGE 7: CONSOLIDATION AND VISUALIZATION\n# =============================================================================\n",
    "generate_qc_report": "\n# =============================================================================\n# STAGE 8: FINAL OUTPUTS AND PACKAGING\n# =============================================================================\n",
}

def extract_rules(content):
    """Extract all rules from the Snakefile content."""
    rules = {}
    
    # Find the preamble (everything before the first rule)
    first_rule_match = re.search(r'^rule\s+\w+:', content, re.MULTILINE)
    if first_rule_match:
        preamble = content[:first_rule_match.start()]
        rules['PREAMBLE'] = preamble.rstrip() + "\n"
        content = content[first_rule_match.start():]
    
    # Regular expression to match rules
    # This matches from "rule name:" to the next "rule" at the start of a line or EOF
    rule_pattern = r'^(rule\s+(\w+):.*?)(?=^rule\s+\w+:|\Z)'
    
    matches = re.finditer(rule_pattern, content, re.MULTILINE | re.DOTALL)
    
    for match in matches:
        rule_content = match.group(1)
        rule_name = match.group(2)
        rules[rule_name] = rule_content
    
    return rules

def reorganize_snakefile(input_file, output_file):
    """Reorganize the Snakefile according to the defined order."""
    
    # Read the original file
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Extract all rules
    rules = extract_rules(content)
    
    # Check for missing rules
    found_rules = set(rules.keys())
    expected_rules = set(RULE_ORDER)
    
    missing_in_order = found_rules - expected_rules
    missing_in_file = expected_rules - found_rules - {'PREAMBLE'}
    
    if missing_in_file:
        print(f"Warning: Rules in order list but not found in file: {missing_in_file}")
    
    if missing_in_order:
        print(f"Warning: Rules found in file but not in order list: {missing_in_order}")
        print("These will be added at the end")
    
    # Build the reorganized content
    reorganized = []
    
    for rule_name in RULE_ORDER:
        if rule_name in rules:
            # Add section header if applicable
            if rule_name in SECTION_HEADERS:
                reorganized.append(SECTION_HEADERS[rule_name])
            
            # Add the rule content
            reorganized.append(rules[rule_name])
            reorganized.append("\n")
    
    # Add any rules that weren't in the order list at the end
    if missing_in_order:
        reorganized.append("\n# =============================================================================\n")
        reorganized.append("# ADDITIONAL RULES (not in predefined order)\n")
        reorganized.append("# =============================================================================\n\n")
        for rule_name in sorted(missing_in_order):
            if rule_name != 'PREAMBLE':
                reorganized.append(rules[rule_name])
                reorganized.append("\n")
    
    # Write the reorganized file
    with open(output_file, 'w') as f:
        f.write(''.join(reorganized))
    
    print(f"Reorganized Snakefile written to {output_file}")
    print(f"Total rules processed: {len(rules) - 1}")  # -1 for PREAMBLE

if __name__ == "__main__":
    input_file = "Snakefile"
    output_file = "Snakefile.reorganized"
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    reorganize_snakefile(input_file, output_file)
    print("\nTo use the reorganized file:")
    print(f"  mv Snakefile Snakefile.backup")
    print(f"  mv {output_file} Snakefile")
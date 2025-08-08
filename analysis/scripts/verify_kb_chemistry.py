#!/usr/bin/env python3
"""
Verify that kb_chemistry strings correctly parse barcode positions from amplicon sequences.
"""

def parse_kb_chemistry(kb_chem_str, amplicon_seq):
    """Parse kb_chemistry string and extract barcodes from amplicon."""
    # Parse the kb_chemistry string
    # Format: umi_start,umi_end,bc3_end,bc2_start,bc2_end,bc2_end,bc1_start,bc1_end,bc1_end:...
    parts = kb_chem_str.split(':')
    positions = [int(x) for x in parts[0].split(',')]
    
    # The format appears to be:
    # positions[0] = 1 (always starts at 1)
    # positions[1] = UMI end position (10)
    # positions[2] = BC3 end position (18)
    # positions[3] = 1 (always 1)
    # positions[4] = BC2 start position (48 means position 49 in 1-indexed)
    # positions[5] = BC2 end position (56 means position 56 in 1-indexed)
    # positions[6] = 1 (always 1)
    # positions[7] = BC1 start position (78 means position 79 in 1-indexed)
    # positions[8] = BC1 end position (86 means position 86 in 1-indexed)
    
    umi_start = positions[0] - 1  # Convert to 0-indexed
    umi_end = positions[1]
    bc3_start = umi_end
    bc3_end = positions[2]
    bc2_start = positions[4]  # Already 0-indexed in the string!
    bc2_end = positions[5]  # Inclusive end position
    bc1_start = positions[7]  # Already 0-indexed in the string!
    bc1_end = positions[8]  # Inclusive end position
    
    # Calculate actual positions
    umi_len = umi_end - umi_start
    bc3_start = umi_end  # BC3 starts right after UMI
    bc3_len = bc3_end - bc3_start
    bc2_len = bc2_end - bc2_start
    bc1_len = bc1_end - bc1_start
    
    # Extract sequences
    umi = amplicon_seq[umi_start:umi_end]
    bc3 = amplicon_seq[bc3_start:bc3_end]
    bc2 = amplicon_seq[bc2_start:bc2_end]
    bc1 = amplicon_seq[bc1_start:bc1_end]
    
    return {
        'UMI': (umi_start, umi_end, umi),
        'BC3': (bc3_start, bc3_end, bc3),
        'BC2': (bc2_start, bc2_end, bc2),
        'BC1': (bc1_start, bc1_end, bc1)
    }

def verify_chemistry(name, amplicon, kb_chemistry, expected_patterns):
    """Verify that extracted barcodes match expected patterns."""
    print(f"\n{name}:")
    print(f"Amplicon: {amplicon}")
    print(f"kb_chemistry: {kb_chemistry}")
    
    result = parse_kb_chemistry(kb_chemistry, amplicon)
    
    all_correct = True
    for bc_name, (start, end, extracted) in result.items():
        expected = expected_patterns.get(bc_name, '')
        is_correct = extracted == expected
        all_correct = all_correct and is_correct
        
        status = "✓" if is_correct else "✗"
        print(f"  {bc_name}: positions {start+1:2d}-{end:2d} = '{extracted}' {status}")
        if not is_correct and expected:
            print(f"       Expected: '{expected}'")
    
    return all_correct

# Define amplicon sequences and chemistry strings
v2_amplicon = "NNNNNNNNNN33333333GTGGCCGATGTTTCGCATCGGCGTACGACT22222222ATCCACGTGCTTGAGACTGTGG11111111"
v2_chemistry = "1,10,18,1,48,56,1,78,86:1,0,10:0,0,0"
v2_expected = {
    'UMI': 'NNNNNNNNNN',
    'BC3': '33333333',
    'BC2': '22222222',
    'BC1': '11111111'
}

v3_amplicon = "NNNNNNNNNN33333333ATGAGGGGTCAG22222222TCCAACCACCTC11111111"
v3_chemistry = "1,10,18,1,30,38,1,50,58:1,0,10:0,0,0"
v3_expected = {
    'UMI': 'NNNNNNNNNN',
    'BC3': '33333333',
    'BC2': '22222222',
    'BC1': '11111111'
}

# Verify both chemistries
v2_correct = verify_chemistry("Chemistry v2", v2_amplicon, v2_chemistry, v2_expected)
v3_correct = verify_chemistry("Chemistry v3", v3_amplicon, v3_chemistry, v3_expected)

print("\n" + "="*50)
if v2_correct and v3_correct:
    print("SUCCESS: All barcode positions correctly parsed!")
else:
    print("ERROR: Some barcode positions are incorrect!")
    exit(1)
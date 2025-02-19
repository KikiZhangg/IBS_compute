#!/usr/bin/env python3

import sys

def parse_23andme(file_path):
    snp_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            # Skip comment or blank lines
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split()
            rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
            
            # Keep only autosomal chromosomes 
            if chrom.isdigit():
                chrom_num = int(chrom)
                if 1 <= chrom_num <= 22:
                    # Skip partial/no-calls; must be exactly two nucleotides
                    if len(genotype) == 2 and '-' not in genotype and '?' not in genotype:
                        # Store in dict
                        snp_dict[rsid] = genotype.upper()
    return snp_dict


def ibs_score(geno1, geno2):
    """
    Return 0, 1, or 2 depending on how many alleles match between geno1 and geno2.
    """
    sorted_g1 = ''.join(sorted(geno1))
    sorted_g2 = ''.join(sorted(geno2))
    
    score = sum(1 for i in range(2) if sorted_g1[i] == sorted_g2[i])
    return score


def compute_ibs_score(file1, file2):
    # Parse both files
    dict1 = parse_23andme(file1)
    dict2 = parse_23andme(file2)
    
    # Find the set of rsids in both
    common_rsids = set(dict1.keys()) & set(dict2.keys())
    
    total_score = 0
    snp_count = 0
    
    for rsid in common_rsids:
        geno1 = dict1[rsid]
        geno2 = dict2[rsid]
        
        if len(geno1) == 2 and len(geno2) == 2:
            score = ibs_score(geno1, geno2)
            total_score += score
            snp_count += 1
    
    # Max possible per SNP is 2
    ibs_similarity = total_score / (2 * snp_count)
    return ibs_similarity



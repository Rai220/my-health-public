#!/usr/bin/env python3
"""
Kinship analysis: determine genetic relatedness between individuals.

Compares genotypes using IBS (Identity by State) statistics:
- IBS0: opposite homozygotes (AA vs BB) — very rare in siblings, ~0 in parent-child
- IBS1: sharing exactly one allele
- IBS2: identical genotypes

For siblings: IBS0 is much lower than for unrelated individuals.
Expected IBD sharing for siblings: P(IBD0)≈0.25, P(IBD1)≈0.50, P(IBD2)≈0.25
"""

import sys
from collections import defaultdict


def parse_genotype_file(filepath):
    """Parse a genotype file into dict of {(chrom, pos): genotype} and {rsid: genotype}."""
    by_pos = {}
    by_rsid = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            rsid, chrom, pos, geno = parts[0], parts[1], parts[2], parts[3]

            # Skip indels and no-calls
            if geno in ('--', 'DD', 'DI', 'II', '') or '/' in geno or len(geno) != 2:
                continue
            # Skip non-standard alleles
            if not all(c in 'ACGT' for c in geno):
                continue

            # Normalize genotype (sort alleles)
            geno_norm = ''.join(sorted(geno))

            key = (str(chrom), str(pos))
            by_pos[key] = geno_norm
            if rsid != '.':
                by_rsid[rsid] = geno_norm

    return by_pos, by_rsid


def compute_ibs(geno1, geno2):
    """
    Compute IBS state between two genotypes.
    Returns 0, 1, or 2 (number of shared alleles).
    """
    a1, a2 = geno1[0], geno1[1]
    b1, b2 = geno2[0], geno2[1]

    # Count shared alleles
    alleles1 = [a1, a2]
    alleles2 = [b1, b2]

    shared = 0
    used = [False, False]
    for a in alleles1:
        for j, b in enumerate(alleles2):
            if not used[j] and a == b:
                shared += 1
                used[j] = True
                break

    return shared


def compare_pair(data1, data2, label1, label2, by='rsid'):
    """Compare two individuals and compute IBS statistics."""
    if by == 'rsid':
        dict1, dict2 = data1[1], data2[1]  # by_rsid
    else:
        dict1, dict2 = data1[0], data2[0]  # by_pos

    common_keys = set(dict1.keys()) & set(dict2.keys())

    ibs_counts = {0: 0, 1: 0, 2: 0}
    het_both = 0  # both heterozygous
    hom_match = 0  # both same homozygous

    # Track chromosome-level stats
    chr_ibs0 = defaultdict(int)
    chr_total = defaultdict(int)

    for key in common_keys:
        g1 = dict1[key]
        g2 = dict2[key]

        ibs = compute_ibs(g1, g2)
        ibs_counts[ibs] += 1

        if g1[0] == g1[1] and g2[0] == g2[1] and g1 == g2:
            hom_match += 1
        if g1[0] != g1[1] and g2[0] != g2[1]:
            het_both += 1

        # Track by chromosome if using position keys
        if by == 'pos':
            chrom = key[0]
        else:
            chrom = 'all'
        if ibs == 0:
            chr_ibs0[chrom] += 1
        chr_total[chrom] += 1

    total = sum(ibs_counts.values())

    print(f"\n{'='*60}")
    print(f"СРАВНЕНИЕ: {label1} vs {label2}")
    print(f"{'='*60}")
    print(f"Общих SNP для сравнения: {total:,}")
    print(f"\nIBS статистика (Identity by State):")
    print(f"  IBS0 (противоположные гомозиготы): {ibs_counts[0]:>8,} ({100*ibs_counts[0]/total:.2f}%)")
    print(f"  IBS1 (общий 1 аллель):             {ibs_counts[1]:>8,} ({100*ibs_counts[1]/total:.2f}%)")
    print(f"  IBS2 (идентичные генотипы):         {ibs_counts[2]:>8,} ({100*ibs_counts[2]/total:.2f}%)")
    print(f"\nДополнительно:")
    print(f"  Оба гомозиготы одинаковые:          {hom_match:>8,} ({100*hom_match/total:.2f}%)")
    print(f"  Оба гетерозиготы:                   {het_both:>8,} ({100*het_both/total:.2f}%)")

    # Relationship estimation
    ibs0_rate = ibs_counts[0] / total
    ibs2_rate = ibs_counts[2] / total

    print(f"\n--- Оценка родства ---")

    print(f"  IBS0 rate: {ibs0_rate:.4f}")
    if ibs0_rate < 0.001:
        print(f"  → Очень низкий IBS0: возможно клоны/близнецы или родитель-ребёнок")
    elif ibs0_rate < 0.04:
        print(f"  → Низкий IBS0: совместимо с родством 1-й степени (родитель-ребёнок или сиблинги)")
        if ibs2_rate > 0.70:
            print(f"  → Высокий IBS2 ({ibs2_rate:.2%}): ближе к родитель-ребёнок")
        else:
            print(f"  → IBS2 = {ibs2_rate:.2%}: ближе к сиблингам (братья/сёстры)")
    elif ibs0_rate < 0.08:
        print(f"  → Умеренный IBS0: совместимо с родством 2-й степени (дядя-племянник, полусиблинги)")
    else:
        print(f"  → Высокий IBS0: вероятно не близкие родственники")

    return ibs_counts, total


def main():
    """
    Example usage. Replace paths with your actual genotype files.

    Genotype file format: TSV with columns: rsid, chromosome, position, genotype
    Example:
        rs9701055   1   630053   CC
        rs3131972   1   817341   AG
    """
    if len(sys.argv) < 3:
        print("Использование: python kinship_analysis.py <file1.txt> <file2.txt> [label1] [label2]")
        print("  file1.txt, file2.txt — файлы генотипов (TSV: rsid, chrom, pos, genotype)")
        print("  label1, label2 — имена для вывода (по умолчанию: Person1, Person2)")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    label1 = sys.argv[3] if len(sys.argv) > 3 else "Person1"
    label2 = sys.argv[4] if len(sys.argv) > 4 else "Person2"

    print(f"Загрузка генотипов {label1}...")
    data1 = parse_genotype_file(file1)
    print(f"  {label1}: {len(data1[1]):,} SNP с rsid, {len(data1[0]):,} по позиции")

    print(f"Загрузка генотипов {label2}...")
    data2 = parse_genotype_file(file2)
    print(f"  {label2}: {len(data2[1]):,} SNP с rsid, {len(data2[0]):,} по позиции")

    # Compare by rsid
    compare_pair(data1, data2, label1, label2, by='rsid')

    # Compare by position
    print(f"\n\n--- Дополнительно: сравнение по позиции ---")
    compare_pair(data1, data2, label1, label2, by='pos')

    print(f"\nСправочные значения IBS0:")
    print(f"  Однояйцевые близнецы:     ~0%")
    print(f"  Родитель-ребёнок:          ~0%")
    print(f"  Полные сиблинги (братья):  ~2-5%")
    print(f"  Полусиблинги:              ~5-8%")
    print(f"  Неродственные:             ~10-15%")


if __name__ == '__main__':
    main()

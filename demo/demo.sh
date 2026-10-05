#!/bin/bash
# ──────────────────────────────────────────────────────────
# LongPhase demo — chr1:121.5-122.0M
# ──────────────────────────────────────────────────────────
set -euo pipefail

LONGPHASE=../longphase
THREADS=4

step() { echo -e "\n\033[1;36m════ $1 ════\033[0m"; }

step "Step 1: Methylation calling"
$LONGPHASE modcall -t $THREADS \
    -b demo.bam \
    -r demo_ref.fa.gz \
    -s demo_snv.vcf.gz \
    -o testing_modcall

step "Step 2a: Phasing (SNV-only, --disableGNN)"
$LONGPHASE phase --ont --disableGNN -t $THREADS \
    -s demo_snv.vcf.gz \
    -r demo_ref.fa.gz \
    -b demo.bam \
    -o testing_snvOnly

step "Step 2b: Phasing (cophasing, --disableGNN)"
$LONGPHASE phase --indels --ont --disableGNN -t $THREADS \
    -s demo_snv.vcf.gz \
    -r demo_ref.fa.gz \
    -b demo.bam \
    --sv-file demo_sv.vcf.gz \
    --mod-file testing_modcall.vcf \
    -o testing_cophasing

step "Step 3a: Phasing (SNV-only, with GNN, the default)"
$LONGPHASE phase --ont -t $THREADS \
    -s demo_snv.vcf.gz \
    -r demo_ref.fa.gz \
    -b demo.bam \
    -o testing_snvOnly_gnn

step "Step 3b: Phasing (cophasing, with GNN, the default)"
$LONGPHASE phase --indels --ont -t $THREADS \
    -s demo_snv.vcf.gz \
    -r demo_ref.fa.gz \
    -b demo.bam \
    --sv-file demo_sv.vcf.gz \
    --mod-file testing_modcall.vcf \
    -o testing_cophasing_gnn

step "Step 4: Evaluate"
for LABEL in testing_snvOnly testing_snvOnly_gnn testing_cophasing testing_cophasing_gnn; do
    $LONGPHASE compare -t $THREADS --ignore-sample-name \
        demo_bench.vcf.gz ${LABEL}.vcf \
        -o sw_${LABEL} 2>&1 | tail -1
done

step "Results"
printf "\n  %-30s %7s %10s %10s\n" "VCF" "SNV_SW" "Phased_SNV" "Hamming%"
printf "  %-30s %7s %10s %10s\n" "──────────────────────────────" "───────" "──────────" "──────────"
for LABEL in testing_snvOnly testing_snvOnly_gnn testing_cophasing testing_cophasing_gnn; do
    TSV="sw_${LABEL}.tsv"
    [[ -f "$TSV" ]] && awk -F'\t' '/^###/ && !/###Sample/ {
        printf "  %-30s %7d %10d %9.3f%%\n", "'"$LABEL"'", $6, $2, $8
    }' "$TSV"
done
echo ""

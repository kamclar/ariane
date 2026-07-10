# Priority Queue Synthesis

Generated: 2026-06-23T12:42:38

## Purpose

This analysis operationalizes VUS prioritization as overlapping manual-review queues. The queues do not add ACMG/ENIGMA evidence and do not reclassify variants. They only label why a VUS may be worth opening earlier.

## Queue Counts

| Queue | VUS count |
| --- | ---: |
| neighborhood context | 6649 |
| low information | 5379 |
| splice/RNA | 3291 |
| functional | 836 |
| near threshold | 637 |
| public assertion | 72 |
| evidence conflict | 30 |

## Number Of Queues Per VUS

| Queue count | VUS count |
| ---: | ---: |
| 1 | 1385 |
| 2 | 4985 |
| 3 | 1603 |
| 4 | 176 |
| 5 | 4 |
| 6 | 1 |

## Interpretation

A high-priority VUS can be high priority for different reasons. Some variants are near an ACMG/ENIGMA threshold, some are splice/RNA candidates, some are useful because public assertions disagree with the local Module 1 output, and many are low-information variants where absence from population data is the main signal.

The most useful manual-review queue is therefore not a universal top-N list. It depends on the curation question: resolving conflicts, finding RNA candidates, reviewing functional assay metadata, or selecting variants close to an accepted classification threshold.

## Output Files

- `tables/priority_queue_synthesis/priority_queue_synthesis_variants.csv`
- `tables/priority_queue_synthesis/priority_queue_counts.csv`
- `tables/priority_queue_synthesis/priority_queue_overlap.csv`
- `tables/priority_queue_synthesis/priority_queue_by_action_group.csv`
- `tables/priority_queue_synthesis/priority_queue_top_multiqueue_examples.csv`
- `plots/30_priority_queue_synthesis/priority_queue_counts.svg`
- `plots/30_priority_queue_synthesis/priority_queue_overlap.svg`

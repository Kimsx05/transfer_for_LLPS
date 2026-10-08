# CellChat consensus sensitivity

- PRIMARY is unchanged: raw sample net$weight, median including eligible real zeros, n_eligible>=3 and prevalence>=0.30.
- For each sample with positive eligible network total, normalized_weight=weight/sum(all eligible sample weights). Eligibility still requires source and target >=30 cells.
- Mathematically undefined zero-total samples: GSM4006647_BC4. They supply NA and are omitted from the NORMALIZED effective denominator only; their real raw zeros remain in the primary network. No normalized zero is invented.
- Normalizable samples=58/59. n_structurally_eligible, n_undefined_total_samples and normalized n_eligible are all retained per edge.
- Spearman edge ranking correlation=0.8567.
- Raw positive edges=182; normalized=182; topology mismatches=0.
- Exact same ordered top-3 external targets=1/7.

## Biological interpretation

- The broad architecture is stable, while relative fibro magnitudes and target ordering are scale-sensitive. C5 has the largest raw total outgoing strength (52.86); sample-normalized cohort edge totals place C0 first (0.566), C1 second (0.375), and C5 third (0.293). C5 also shifts from first to third for total incoming strength.
- The strongest external raw edge is C0→Epithelial (4.146; 16 eligible samples, detected in all 16). C5→Epithelial and C5→Endothelial are also strong and recurrent. Raw versus normalized top receiver switches include C1 (Endothelial→Epithelial) and C5 (Epithelial→Endothelial). Epithelial and Endothelial remain central fibro receiver relationships.
- This is a material qualification of subtype dominance: do not claim C5 is unconditionally the most active fibro subtype from raw row totals alone. Raw weight measures absolute inferred probability sums; normalized weight measures relative communication within each sample. Different subtype eligible-sample sets and sample total signal scales contribute to the change.
- Row/column totals sum edge medians calculated from different eligible sample subsets. They are descriptive summaries, not a cohort-wide probability mass, and not cell-number-adjusted activity. C6 retains its low-coverage flag (3 eligible samples).
- Primary raw network and cutoffs remain unchanged. NicheNet follow-up should retain C0 and C5 and both epithelial/endothelial receiver programs in the design queue; no gene program or ligand activity is selected here.

## External target ordering

Key: <source>
     source                     target_raw              target_normalized
     <char>                         <char>                         <char>
1: Fibro_C0 Epithelial;Endothelial;Myeloid        Epithelial;Myeloid;T_NK
2: Fibro_C1    Endothelial;Epithelial;Mast Epithelial;Endothelial;Myeloid
3: Fibro_C2    Epithelial;Endothelial;Mast Epithelial;Myeloid;Endothelial
4: Fibro_C3    Epithelial;Endothelial;Mast    Epithelial;Endothelial;T_NK
5: Fibro_C4    Epithelial;Endothelial;Mast    Epithelial;Endothelial;T_NK
6: Fibro_C5    Epithelial;Endothelial;T_NK Endothelial;Epithelial;Myeloid
7: Fibro_C6       Endothelial;T_NK;Myeloid       Endothelial;T_NK;Myeloid
   same_top3_targets
              <lgcl>
1:             FALSE
2:             FALSE
3:             FALSE
4:             FALSE
5:             FALSE
6:             FALSE
7:              TRUE

## Relative C0–C6 patterns

   fibro_subtype outgoing_raw incoming_raw outgoing_normalized
          <char>        <num>        <num>               <num>
1:      Fibro_C0     36.67198    19.086497           0.5655783
2:      Fibro_C1     41.67282    31.832459           0.3749236
3:      Fibro_C2     28.92338     9.670982           0.2098982
4:      Fibro_C3     24.55298    20.976765           0.1651374
5:      Fibro_C4     31.22870    25.760898           0.2515140
6:      Fibro_C5     52.85705    36.451029           0.2927409
7:      Fibro_C6     28.22579    21.145242           0.1334268
   incoming_normalized outgoing_to_major_raw outgoing_to_major_normalized
                 <num>                 <num>                        <num>
1:          0.22109531              16.55582                   0.37692193
2:          0.24637053              15.61444                   0.19478260
3:          0.06507293              13.26548                   0.09959206
4:          0.14569734              12.27530                   0.08321014
5:          0.16841143              14.64503                   0.12702978
6:          0.20382970              20.22988                   0.11742555
7:          0.10109105              11.68039                   0.05582569
   incoming_from_major_raw incoming_from_major_normalized outgoing_rank_raw
                     <num>                          <num>             <int>
1:                3.087724                     0.04844375                 3
2:                3.430024                     0.04469430                 2
3:                1.758174                     0.01327075                 5
4:                3.508180                     0.02795636                 7
5:                3.992276                     0.02939354                 4
6:                4.993187                     0.02735100                 1
7:                4.287962                     0.02202707                 6
   outgoing_rank_normalized incoming_rank_raw incoming_rank_normalized
                      <int>             <int>                    <int>
1:                        1                 6                        2
2:                        2                 2                        1
3:                        5                 7                        7
4:                        6                 5                        5
5:                        4                 3                        4
6:                        3                 1                        3
7:                        7                 4                        6

- The initial stricter handling of undefined normalized observations is archived under 00_audit/sensitivity_initial_undefined_total_handling. The final sensitivity omits only samples whose normalization denominator is undefined; the primary network was never recomputed or overwritten.

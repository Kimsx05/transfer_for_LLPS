# CellChat consensus sensitivity

- Primary remains unnormalized per-sample net$weight -> median over all eligible samples including zeros -> n_eligible>=3 and prevalence>=0.30.
- Sensitivity divides each sample weight by the sum of all eligible weights in that same sample; no pooling.
- Samples with undefined/zero network total: 1. Such samples supply NA, not invented normalized zeros.
- Spearman edge-strength ranking correlation: 0.8643.
- Raw positive displayed edges: 182; normalized: 181; topology mismatches: 1.
- Same fibro top-3 external targets: 1/7.
- Retained-by-rule edges with median=0 remain in the audit but cannot have visible positive edge width. This follows the requested zero-inclusive median definition.

## Fibro target comparison

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

## Relative C0–C6 pattern

   fibro_subtype outgoing_raw incoming_raw outgoing_normalized
          <char>        <num>        <num>               <num>
1:      Fibro_C0     36.67198    19.086497           0.5655783
2:      Fibro_C1     41.67282    31.832459           0.3749236
3:      Fibro_C2     28.92338     9.670982           0.2098982
4:      Fibro_C3     24.55298    20.976765           0.1651374
5:      Fibro_C4     31.22870    25.760898           0.2515140
6:      Fibro_C5     52.85705    36.451029           0.2927409
7:      Fibro_C6     28.22579    21.145242           0.1334268
   incoming_normalized outgoing_rank_raw outgoing_rank_normalized
                 <num>             <int>                    <int>
1:          0.22109531                 3                        1
2:          0.24637053                 2                        2
3:          0.06507293                 5                        5
4:          0.14569734                 7                        6
5:          0.16841143                 4                        4
6:          0.20382970                 1                        3
7:          0.10109105                 6                        7
   incoming_rank_raw incoming_rank_normalized
               <int>                    <int>
1:                 6                        2
2:                 2                        1
3:                 7                        7
4:                 5                        5
5:                 3                        4
6:                 1                        3
7:                 4                        6

- Rank/top-target discrepancies are descriptive checks, not statistical significance tests. Incomplete edges caused by failed inference are explicitly unavailable; no failure is converted to a zero.

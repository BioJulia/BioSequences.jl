@test molecular_weight(aa"PKLEQ") == 613.68
@test molecular_weight(aa"ETIWS*") == 634.65

@test molecular_weight(dna"GCAGCCATAG") == 3036.930000000001
@test molecular_weight(dna"TCCCAGACTG", :phosphate, :double) == 6215.960000000001

@test molecular_weight(rna"GGGCGGACCU") == 3214.9
@test molecular_weight(rna"CGAUUUUCGG", :triphosphate, :double) == 6784.700000000001

@test _molecular_weight(dna"GTTGCCCGGC",DNA_WEIGHTS) == 3019.9000000000005
@test _molecular_weight(rna"GUCUGACGCG",RNA_WEIGHTS, :triphosphate) == 3415.8600000000006


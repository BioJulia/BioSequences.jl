# Sources of the molecular weights of each defined amino acid
# AA_A   # Alanine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5950
# AA_R   # Arginine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6322
# AA_N   # Asparagine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6267
# AA_D   # Aspartate                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5960
# AA_C   # Cysteine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5862
# AA_Q   # Glutamine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5961
# AA_E   # Glutamate                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/33032
# AA_G   # Glycine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/750
# AA_H   # Histidine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6274
# AA_I   # Isoleucine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6306
# AA_L   # Leucine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6106
# AA_K   # Lysine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5962
# AA_M   # Methionine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6137
# AA_F   # Phenylalanine               ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6140 
# AA_P   # Proline                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/145742
# AA_S   # Serine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5951
# AA_T   # Threonine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6288
# AA_W   # Tryptophan                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6305
# AA_Y   # Tyrosine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6057
# AA_V   # Valine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6287
# AA_O   # Pyrrolysine                 ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/21873141
# AA_U   # Selenocysteine              ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/25076

# Source of the used water molecular weight = 18.02 gmol for calculations
# Water                                ->  Source = https://www.cmu.edu/gelfand/k12-educational-resources/polymers/what-is-polymer/molecular-weight-calculation.html 

# Creating the array AA_WEIGHTS of length 28 to list all weights. In ambigous cases, let's list -1.0 
# In the function, if the weight is -1, then it will trow an error
# Note that for gap and terminal codon, the value is 18.02 intead of 0 to account for the substraction of water in the function
const AA_WEIGHTS = [89.09, 174.20, 132.12, 133.10, 121.16, 146.14, 147.13, 75.07, 155.15, 131.17, 131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.22, 181.19, 117.15, 255.31, 168.06, -1.0, 131.17, -1.0, -1.0, 18.02, 18.02]

# Defining the function molecular_weight, which accepts an amino acid sequence and return molecular weight in  g/mol. Note that the function also account for lost of water for each peptide bond formation 
# Results are consistent with the one of of this website: https://www.bioinformatics.org/sms/prot_mw.html

"""
    molecular_weight(aa_seq::AASeq) -> Float64

Calculate molecular weight of `aa_seq`in g/mol, i.e. 
Sum of the weights of all the amino acids of the sequence minus water lost per peptide bond. 
Values of each amino acid were obtained from PubChem 2.1, values are listed in the vector{Float64} AA_WEIGHTS.

# Examples
```jldoctest
julia> molecular_weight(aa"PKLEQ")
613.68

julia> molecular_weight(aa"ETIWS*")
634.65
```
"""
function molecular_weight(aa_seq::AASeq)
    weight = 0.0
    for aa in aa_seq
        if AA_WEIGHTS[reinterpret(UInt8, aa) + 1] == -1.0
            throw(ArgumentError("amino acid $aa weight is ambiguous"))
        else 
            aa_weight = AA_WEIGHTS[reinterpret(UInt8, aa) + 1]
            weight += aa_weight
        end
    end
    return weight - ((length(aa_seq) - 1) * 18.02)
end

# Sources of the molecular weights of each defined nucleotide
# DNA Adenine                 dAMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/12599
# RNA Adenine                 AMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6083
# DNA Cytosine                dCMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/13945
# RNA Cytosine                CMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6131  
# DNA Guanine                 dGMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/135398597
# RNA Guanine                 GMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/135398631
# DNA Thymine                 dTMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/9700
# RNA Uracil                  UMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6030
# Calculations following the following guidelines: https://ymc.eu/files/imported/publications/556/documents/YMC-Expert-Tip---How-to-calculate-the-MW-of-nucleic-acids.pdf

# Creating the arrays DNA_WEIGHTS and RNA_WEIGHTS to list all weights. In ambigous cases, let's list -1.0, so if called it trows an error
# The listed values does not account for the loss of water and it will be taken into consideration in the _molecular weight function. The used value for water weight is the same as above
const DNA_WEIGHTS = [0, 331.22, 307.20, -1.0, 347.22, -1.0, -1.0, -1.0, 322.21, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

const RNA_WEIGHTS = [0, 347.22, 323.20, -1.0, 363.22, -1.0, -1.0, -1.0, 324.18 , -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

# Creating 3 functions which can process DNA or RNA seqeunces, with variable five_terminal_state and strand_number
# If not specified the assumption are single stranded and hydroxyl as terminal 5' functional group
# Results are consistent with the one of of this website: https://molbiotools.com/dnacalculator.php

"""
    molecular_weight(rna_seq::LongSequence{DNAAlphabet{4}}, five_terminal_state::Symbol, strand_number::Symbol) -> Float64

Calculate molecular weight of `dna_seq`in g/mol, i.e. Sum of the weights of all the nucleotides of the sequence, 
minus water lost per bond formation. 
Values of nucleotides were obtained from PubChem 2.1 and are listed in the vector{Float64} DNA_WEIGHTS.
Also possible to define the 5' state as (:hydroxyl, :phosphate or :triphosphate) and strand number as (:single or :double).
If not defined, the default values are :hydroxyl and :single

# Examples
```jldoctest
julia> molecular_weight(dna"GCAGCCATAG")
3036.930000000001

julia> mmolecular_weight(dna"TCCCAGACTG", :phosphate, :double)
6215.960000000001
```
"""

function molecular_weight(dna_seq::LongSequence{DNAAlphabet{4}}, five_terminal_state=:hydroxyl, strand_number=:single)
    if strand_number == :single
        return _molecular_weight(dna_seq, DNA_WEIGHTS, five_terminal_state)
    elseif strand_number == :double
        com_dna_seq = complement(dna_seq)
        return _molecular_weight(dna_seq, DNA_WEIGHTS, five_terminal_state) + _molecular_weight(com_dna_seq, DNA_WEIGHTS, five_terminal_state)     
    else
        throw(ArgumentError("Unknown strand_number $strand_number. Must be :single or :double"))
    end
end

"""
    molecular_weight(rna_seq::LongSequence{RNAAlphabet{4}}, five_terminal_state::Symbol, strand_number::Symbol) -> Float64

Calculate molecular weight of `rna_seq`in g/mol, i.e. Sum of the weights of all the nucleotides of the sequence, 
minus water lost per bond formation. 
Values of nucleotides were obtained from PubChem 2.1 and are listed in the vector{Float64} RNA_WEIGHTS.
Also possible to define the 5' state as (:hydroxyl, :phosphate or :triphosphate) and strand number as (:single or :double).
If not defined, the default values are :hydroxyl and :single

# Examples
```jldoctest
julia> molecular_weight(rna"GGGCGGACCU")
3214.9

julia> molecular_weight(rna"CGAUUUUCGG", :triphosphate, :double)
6784.700000000001
```
"""

function molecular_weight(rna_seq::LongSequence{RNAAlphabet{4}}, five_terminal_state=:hydroxyl, strand_number=:single)
    if strand_number == :single
        return _molecular_weight(rna_seq, RNA_WEIGHTS, five_terminal_state)
    elseif strand_number == :double
        com_rna_seq = complement(rna_seq)
        return _molecular_weight(rna_seq, RNA_WEIGHTS, five_terminal_state) + _molecular_weight(com_rna_seq, RNA_WEIGHTS, five_terminal_state)     
    else
        throw(ArgumentError("Unknown strand_number $strand_number. Must be :single or :double"))
    end
end

"""
    _molecular_weight(nucseq::NucSeq, array::Vector{Float64}, five_terminal_state::Symbol) -> Float64

Calculate molecular weight of `rnucseq`in g/mol, i.e. Sum of the weights of all the nucleotides of the sequence, 
minus water lost per bond. Values of nucleotides were obtained from PubChem 2.1. 
Values are listed in the vectors{Float64} RNA_WEIGHTS and DNA_WEIGHTS.
Also possible to define the 5' state as (:hydroxyl, :phosphate or :triphosphate), If not defined, 
hydroxyl is the default option. Not possible to calculate weights for double stranded sequences.

# Examples
```jldoctest
julia> _molecular_weight(dna"GTTGCCCGGC",DNA_WEIGHTS)
3019.9000000000005

julia> _molecular_weight(rna"GUCUGACGCG",RNA_WEIGHTS, :triphosphate)
3415.8600000000006
```
"""

function _molecular_weight(nucseq::NucSeq, array::Vector{Float64}, five_terminal_state=:hydroxyl)
    weight = 0.0
    for nucleotide in nucseq
        if array[reinterpret(UInt8, nucleotide) + 1] == -1.0
            throw(ArgumentError("nucleotide $nucleotide weight is ambiguous"))
        else
            dna_weight = array[reinterpret(UInt8, nucleotide) + 1]
            weight += dna_weight
        end
    end
    if five_terminal_state == :hydroxyl
        weight =  weight - 62
    elseif five_terminal_state == :phosphate
        weight = weight + 18.06
    elseif five_terminal_state == :triphosphate 
        weight = weight + 178
    else
        throw(ArgumentError("Unknown five_terminal_state $five_terminal_state. Must be :hydroxyl, :phosphate or :triphosphate"))  
    end
    return  weight - (length(nucseq) * 18.02)
end
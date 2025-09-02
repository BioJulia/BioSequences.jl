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

# Creating the array AA_WEIGHTS of length 28 to list all weights. In ambigous cases, let's list the string "Undefined", so if called, it returns an error
AA_WEIGHTS = [89.09, 174.20, 132.12, 133.10, 121.16, 146.14, 147.13, 75.07, 155.15, 131.17, 131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.22, 181.19, 117.15, 255.31, 168.06, "Undefined", 131.17, "Undefined", "Undefined", 0, 0]

# Defining the function molecular_weight, which accepts an amino acid sequence and return molecular weight in  g/mol. Note that the function also account for lost of water for each peptide bond formation 
function molecular_weight(aa_seq::AASeq)
    weight = 0
    for aa in aa_seq
        aa_weight = AA_WEIGHTS[reinterpret(UInt8, aa) + 1]
        weight += aa_weight
    end
    return round(weight - ((length(aa_seq) - 1) * 18.02); digits = 3)
end

# Sources of the molecular weights of each defined nucleotide
# DNA Adenine                 dAMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/12599
# RNA Adenine                 AMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6083
# DNA Cytosine                dCMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/13945
# RNA Cytosine                CMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6131  
# DNA Guanine                 dGMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/1135398597
# RNA Guanine                 GMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/135398631
# DNA Thymine                 dTMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/9700
# Calculations following the following guidelines: https://ymc.eu/files/imported/publications/556/documents/YMC-Expert-Tip---How-to-calculate-the-MW-of-nucleic-acids.pdf

# Creating the arrays DNA_WEIGHTS and RNA_WEIGHTS to list all weights. In ambigous cases, let's list the string "Undefined", so if called, it returns an error
# The values already account for loss of water = - 18.02
DNA_WEIGHTS = [0, 313.2, 289.18, "Undefined", 329.2, "Undefined", "Undefined", "Undefined", 304.19, "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined"]

RNA_WEIGHTS = [0, 329.2, 305.18, "Undefined", 345.2, "Undefined", "Undefined", "Undefined", 306.16, "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined", "Undefined"]

# Defining a function that accepts nucleotide sequences (either DNA or RNA) and calculate the molecular weight in g/mol
# Note the additional keywords alphabet, five_terminal_state and strand_number, which have default cases, if not specified
# In case of using double strand, the five_terminal_state of the complmentary strand is assumed to be the same as the input strand
# Function only accepts double as strand_number for DNA sequences
function molecular_weight(nucseq::NucSeq, alphabet=:DNA, five_terminal_state=:hydroxyl, strand_number=:single)
    weight = 0
    if alphabet == :DNA
      for dnucleotide in nucseq
        dna_weight = DNA_WEIGHTS[reinterpret(UInt8, dnucleotide) + 1]
        weight += dna_weight
        end
    elseif alphabet == :RNA
      for nucleotide in nucseq
        rna_weight = RNA_WEIGHTS[reinterpret(UInt8, nucleotide) + 1]
        weight += rna_weight
        end 
    else 
        throw(ArgumentError("Unknown Alphabet. Must be :DNA or :RNA"))
    end
    if five_terminal_state == :hydroxyl
        weight =  weight - 62
    elseif five_terminal_state == :phosphate
        weight = weight + 79
    elseif five_terminal_state == :triphosphate 
        weight = weight + 178
    else
        throw(ArgumentError("Unknown five_terminal_state. Must be :hydroxyl, :phosphate or :triphosphate"))  
    end
    if strand_number == :single
        return round(weight; digits = 3)
    elseif strand_number == :double && alphabet == :DNA
        com_seq = complement(nucseq)
        com_weight = 0
        for c_dnucleotide in com_seq
            com_dna_weight = DNA_WEIGHTS[reinterpret(UInt8, c_dnucleotide) + 1]
            com_weight += com_dna_weight
        end
        if five_terminal_state == :hydroxyl
        com_weight =  com_weight - 62
        elseif five_terminal_state == :phosphate
        com_weight = com_weight + 79
        elseif five_terminal_state == :triphosphate
        com_weight = com_weight + 178 
        end
        return round(com_weight + weight; digits = 3)
    else 
        throw(ArgumentError("Unknown strand_number. Must be :single or :double. Or trying to use double as strand_number with a RNA sequence"))
    end    
end
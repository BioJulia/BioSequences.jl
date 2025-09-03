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
const AA_WEIGHTS = [89.09, 174.20, 132.12, 133.10, 121.16, 146.14, 147.13, 75.07, 155.15, 131.17, 131.17, 146.19, 149.21, 165.19, 115.13, 105.09, 119.12, 204.22, 181.19, 117.15, 255.31, 168.06, -1.0, 131.17, -1.0, -1.0, 0, 0]

# Defining the function molecular_weight, which accepts an amino acid sequence and return molecular weight in  g/mol. Note that the function also account for lost of water for each peptide bond formation 
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

# Creating the arrays DNA_WEIGHTS and RNA_WEIGHTS to list all weights. In ambigous cases, let's list -1.0, so if called
# The listed values does not account for the loss of water and it will be taken into consideration later
# In the function, if the weight is -1, then it will trow an error
const DNA_WEIGHTS = [0, 331.22, 307.20, -1.0, 347.22, -1.0, -1.0, -1.0, 322.21, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

const RNA_WEIGHTS = [0, 347.22, 323.20, -1.0, 363.22, -1.0, -1.0, -1.0, 324.18 , -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0]

# Defining a function that accepts nucleotide sequences (either DNA or RNA) and calculate the molecular weight in g/mol
# Note the additional keywords alphabet, five_terminal_state and strand_number, which have default cases, if not specified
# In case of using double strand, the five_terminal_state of the complmentary strand is assumed to be the same as the input strand
# Function only accepts double as strand_number for DNA sequences
function molecular_weight(nucseq::NucSeq, alphabet=:DNA, five_terminal_state=:hydroxyl, strand_number=:single)
    weight = 0.0
    if alphabet == :DNA
      for dnucleotide in nucseq
        if DNA_WEIGHTS[reinterpret(UInt8, dnucleotide) + 1] == -1.0
            throw(ArgumentError("nucleotide $dnucleotide weight is ambiguous"))
        else
            dna_weight = DNA_WEIGHTS[reinterpret(UInt8, dnucleotide) + 1]
            weight += dna_weight
        end
       end
    elseif alphabet == :RNA
      for nucleotide in nucseq
        if RNA_WEIGHTS[reinterpret(UInt8, nucleotide) + 1] == -1.0
            throw(ArgumentError("nucleotide $nucleotide weight is ambiguous"))
        else
            rna_weight = RNA_WEIGHTS[reinterpret(UInt8, nucleotide) + 1]
            weight += rna_weight
        end
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
        return weight - (length(nucseq) * 18.02) 
    elseif strand_number == :double && alphabet == :DNA
        com_seq = complement(nucseq)
        com_weight = 0.0
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
        return com_weight + weight - (2 * length(nucseq) * 18.02)
    else 
        throw(ArgumentError("Unknown strand_number. Must be :single or :double. Or trying to use double as strand_number with a RNA sequence"))
    end    
end
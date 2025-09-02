# Sources of the molecular weights of each defined amino acid in g/mol
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

# Definning const of molecular weight of monophosphate in g/mol to account for the loss of water and PPi during synthesis
const DNA_A_MW = 331.22 - WATER_MW          # DNA Adenine                 dAMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/12599
const RNA_A_MW = 347.22 - WATER_MW          # RNA Adenine                 AMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6083
const DNA_C_MW = 307.20 - WATER_MW          # DNA Cytosine                dCMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/13945
const RNA_C_MW = 323.20 - WATER_MW          # RNA Cytosine                CMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6131  
const DNA_G_MW = 347.22 - WATER_MW          # DNA Guanine                 dGMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/1135398597
const RNA_G_MW = 363.22 - WATER_MW          # RNA Guanine                 GMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/135398631
const DNA_T_MW = 322.21 - WATER_MW          # DNA Thymine                 dTMP ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/9700
const RNA_U_MW = 324.18 - WATER_MW          # RNA Uracil                  UMP  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6030

# Calculations following the following guidelines: https://ymc.eu/files/imported/publications/556/documents/YMC-Expert-Tip---How-to-calculate-the-MW-of-nucleic-acids.pdf
# Defining a function, which accepts DNA sequences
function molecular_weight(aa_seq::NucSeq{4, DNAAlphabet{4}}, unit=:"gmol", five_terminal_state=:"Hydroxyl")
    AA_vector = Vector{Float64}(undef, 0)
    i = 1
    while i <= length(aa_seq)
        if aa_seq[i] == DNA_A
            push!(AA_vector, DNA_A_MW)
        elseif aa_seq[i] == DNA_C
            push!(AA_vector, DNA_C_MW)
        elseif aa_seq[i] == DNA_G
            push!(AA_vector, DNA_G_MW)
        elseif aa_seq[i] == DNA_T
            push!(AA_vector, DNA_T_MW)                       
        else
            println("Error due to unknown nucleotide at $i position")
            return nothing
        end
        i = i+1
    end
        if unit == "gmol" && five_terminal_state == "Phosphate"
            return round(sum(AA_vector) + 79; digits = 3)
        elseif unit == "gmol" && five_terminal_state == "Hydroxyl"
            return round(sum(AA_vector) - 62; digits = 3)
        elseif unit == "gmol" && five_terminal_state == "Triphosphate"
            return round(sum(AA_vector) + 178; digits = 3)
        elseif unit == "kDa" && five_terminal_state == "Phosphate"
            return (round(sum(AA_vector) + 79; digits = 3))/1000
        elseif unit == "kDa" && five_terminal_state == "Hydroxyl"
            return (round(sum(AA_vector) - 62; digits = 3))/1000
        elseif unit == "kDa" && five_terminal_state == "Triphosphate"
            return (round(sum(AA_vector) + 178; digits = 3))/1000        
        else 
            println("Error")
            return nothing
        end
end

# Defining a function, which accepts a RNA sequences
function molecular_weight(aa_seq::NucSeq{4, RNAAlphabet{4}}, unit=:"gmol", five_terminal_state=:"Triphosphate")
    AA_vector = Vector{Float64}(undef, 0)
    i = 1
    while i <= length(aa_seq)
        if aa_seq[i] == RNA_A
            push!(AA_vector, RNA_A_MW)
        elseif aa_seq[i] == RNA_C
            push!(AA_vector, RNA_C_MW)
        elseif aa_seq[i] == RNA_G
            push!(AA_vector, RNA_G_MW)
        elseif aa_seq[i] == RNA_U
            push!(AA_vector, RNA_U_MW)                       
        else
            println("Error due to unknown nucleotide at $i position")
            return nothing
        end
        i = i+1
    end
        if unit == "gmol" && five_terminal_state == "Triphosphate"
            return round(sum(AA_vector) + 178; digits = 3)
        elseif unit == "gmol" && five_terminal_state == "Phosphate"
            return round(sum(AA_vector) + 79; digits = 3)            
        elseif unit == "gmol" && five_terminal_state == "Hydroxyl"
            return round(sum(AA_vector) - 62; digits = 3)
        elseif unit == "kDa" && five_terminal_state == "Triphosphate"
            return (round(sum(AA_vector) - 62; digits = 3))/1000  
        elseif unit == "kDa" && five_terminal_state == "Phosphate"
            return (round(sum(AA_vector) + 79; digits = 3))/1000
        elseif unit == "kDa" && five_terminal_state == "Hydroxyl"
            return (round(sum(AA_vector) - 62; digits = 3))/1000        
        else 
            println("Error")
            return nothing
        end
end
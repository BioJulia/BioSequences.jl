# Definning const of molecular weight of canonical amino acids in g/mol
const AA_A_MW = 89.09                     # Alanine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5950
const AA_R_MW = 174.20                    # Arginine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6322
const AA_N_MW = 132.12                    # Asparagine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6267
const AA_D_MW = 133.10                    # Aspartate                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5960
const AA_C_MW = 121.16                    # Cysteine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5862
const AA_Q_MW = 146.14                    # Glutamine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5961
const AA_E_MW = 147.13                    # Glutamate                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/33032
const AA_G_MW = 75.07                     # Glycine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/750
const AA_H_MW = 155.15                    # Histidine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6274
const AA_I_MW = 131.17                    # Isoleucine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6306
const AA_L_MW = 131.17                    # Leucine                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6106
const AA_K_MW = 146.19                    # Lysine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5962
const AA_M_MW = 149.21                    # Methionine                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6137
const AA_F_MW = 165.19                    # Phenylalanine               ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6140 
const AA_P_MW = 115.13                    # Proline                     ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/145742
const AA_S_MW = 105.09                    # Serine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/5951
const AA_T_MW = 119.12                    # Threonine                   ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6288
const AA_W_MW = 204.22                    # Tryptophan                  ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6305
const AA_Y_MW = 181.19                    # Tyrosine                    ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6057
const AA_V_MW = 117.15                    # Valine                      ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/6287
# Defining water molecular weight (required for calculation of water lost)
const WATER_MW = 18.02                    # Water                       ->  Source = https://www.cmu.edu/gelfand/k12-educational-resources/polymers/what-is-polymer/molecular-weight-calculation.html

# Non canonical / not defined amino acids in g/mol
const AA_O_MW_MW  = 255.31                    # Pyrrolysine                 ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/21873141
const AA_U_MW_MW  = 168.06                    # Selenocysteine              ->  Source = https://pubchem.ncbi.nlm.nih.gov/compound/25076
const AA_B_MW_MW  = (AA_N_MW + AA_D_MW)/2     # Aspartate or Asparagine     ->  Average between Aspartic acid and Asparagine
const AA_J_MW_MW  = 131.17                    # Leucine or Isoleucine       ->  Both have the same molecular weight
const AA_Z_MW_MW  = (AA_Q_MW + AA_D_MW)/2     # Glutamine or Glutamine      ->  Average between Glutamine and Glutamine
const AA_X_MW_MW  = 110                       # Any amino acid (average)    ->  Source = http://thermofisher.com/de/de/home/references/ambion-tech-support/rna-tools-and-calculators/proteins-and-amino-acids.html
const AA_O_Gap_MW = 0                         # Consider gap as 0       

# Defining the function molecular_weight, which accepts a AA sequence and units ("gmol", "KDa") as input and return molecular weight in either g/mol or kDa. If not specified it returns in g/mol
function molecular_weight(aa_seq::LongAA, unit=:"gmol")
    AA_vector = Vector{Float64}(undef, 0)
    i = 1
    while i <= length(aa_seq)
        if aa_seq[i] == AA_A
            push!(AA_vector, AA_A_MW)
        elseif aa_seq[i] == AA_R
            push!(AA_vector, AA_R_MW)
        elseif aa_seq[i] == AA_N
            push!(AA_vector, AA_N_MW)
        elseif aa_seq[i] == AA_D
            push!(AA_vector, AA_D_MW)
        elseif aa_seq[i] == AA_C
            push!(AA_vector, AA_C_MW)
        elseif aa_seq[i] == AA_Q
            push!(AA_vector, AA_Q_MW)
        elseif aa_seq[i] == AA_E
            push!(AA_vector, AA_E_MW)
        elseif aa_seq[i] == AA_G
            push!(AA_vector, AA_G_MW)
        elseif aa_seq[i] == AA_H
            push!(AA_vector, AA_H_MW)  
        elseif aa_seq[i] == AA_I
            push!(AA_vector, AA_I_MW)
        elseif aa_seq[i] == AA_L
            push!(AA_vector, AA_L_MW)    
        elseif aa_seq[i] == AA_K
            push!(AA_vector, AA_K_MW)
        elseif aa_seq[i] == AA_M
            push!(AA_vector, AA_M_MW)                                            
        elseif aa_seq[i] == AA_F
            push!(AA_vector, AA_F_MW)
        elseif aa_seq[i] == AA_P
            push!(AA_vector, AA_P_MW)
        elseif aa_seq[i] == AA_S
            push!(AA_vector, AA_S_MW)
        elseif aa_seq[i] == AA_T
            push!(AA_vector, AA_T_MW)   
        elseif aa_seq[i] == AA_W
            push!(AA_vector, AA_W_MW)
        elseif aa_seq[i] == AA_Y
            push!(AA_vector, AA_Y_MW)
        elseif aa_seq[i] == AA_V
            push!(AA_vector, AA_V_MW)
        elseif aa_seq[i] == AA_O
            push!(AA_vector, AA_O_MW_MW)
        elseif aa_seq[i] == AA_U
            push!(AA_vector, AA_U_MW_MW)
        elseif aa_seq[i] == AA_B
            push!(AA_vector, AA_B_MW_MW)
        elseif aa_seq[i] == AA_J
            push!(AA_vector, AA_J_MW_MW)
        elseif aa_seq[i] == AA_Z
            push!(AA_vector, AA_Z_MW_MW)
        elseif aa_seq[i] == AA_X
            push!(AA_vector, AA_X_MW_MW)
        elseif aa_seq[i] == AA_O_Gap
            push!(AA_vector, AA_O_Gap_MW)                           
        else
            println("Error due to unknown amino acid at $i position")
            return nothing
        end
        i = i+1
    end
        if unit == "gmol"
            return round(sum(AA_vector) - ((length(aa_seq) - 1) * WATER_MW); digits = 3)
        elseif unit == "kDa"
            return round(((sum(AA_vector) - ((length(aa_seq) - 1) * WATER_MW))/1000); digits = 3)
        else 
            return nothing
        end
end

# Defining a function, which accepts a AA sequence and print the molecular weight in g/mol and kDa from a given n_start position in the sequence until a n_end given position in the sequence
function molecular_weight(aa_seq::LongAA, n_start::Integer, n_end::Integer)
    AA_vector = Vector{Float64}(undef, 0)
    i = n_start
    while i <= n_end
        if aa_seq[i] == AA_A
            push!(AA_vector, AA_A_MW)
        elseif aa_seq[i] == AA_R
            push!(AA_vector, AA_R_MW)
        elseif aa_seq[i] == AA_N
            push!(AA_vector, AA_N_MW)
        elseif aa_seq[i] == AA_D
            push!(AA_vector, AA_D_MW)
        elseif aa_seq[i] == AA_C
            push!(AA_vector, AA_C_MW)
        elseif aa_seq[i] == AA_Q
            push!(AA_vector, AA_Q_MW)
        elseif aa_seq[i] == AA_E
            push!(AA_vector, AA_E_MW)
        elseif aa_seq[i] == AA_G
            push!(AA_vector, AA_G_MW)
        elseif aa_seq[i] == AA_H
            push!(AA_vector, AA_H_MW)  
        elseif aa_seq[i] == AA_I
            push!(AA_vector, AA_I_MW)
        elseif aa_seq[i] == AA_L
            push!(AA_vector, AA_L_MW)    
        elseif aa_seq[i] == AA_K
            push!(AA_vector, AA_K_MW)
        elseif aa_seq[i] == AA_M
            push!(AA_vector, AA_M_MW)                                            
        elseif aa_seq[i] == AA_F
            push!(AA_vector, AA_F_MW)
        elseif aa_seq[i] == AA_P
            push!(AA_vector, AA_P_MW)
        elseif aa_seq[i] == AA_S
            push!(AA_vector, AA_S_MW)
        elseif aa_seq[i] == AA_T
            push!(AA_vector, AA_T_MW)   
        elseif aa_seq[i] == AA_W
            push!(AA_vector, AA_W_MW)
        elseif aa_seq[i] == AA_Y
            push!(AA_vector, AA_Y_MW)
        elseif aa_seq[i] == AA_V
            push!(AA_vector, AA_V_MW)
        elseif aa_seq[i] == AA_O
            push!(AA_vector, AA_O_MW_MW)
        elseif aa_seq[i] == AA_U
            push!(AA_vector, AA_U_MW_MW)
        elseif aa_seq[i] == AA_B
            push!(AA_vector, AA_B_MW_MW)
        elseif aa_seq[i] == AA_J
            push!(AA_vector, AA_J_MW_MW)
        elseif aa_seq[i] == AA_Z
            push!(AA_vector, AA_Z_MW_MW)
        elseif aa_seq[i] == AA_X
            push!(AA_vector, AA_X_MW_MW)
        elseif aa_seq[i] == AA_O_Gap
            push!(AA_vector, AA_O_Gap_MW)                           
        else
            println("Error due to unknown amino acid at $i position")
            return nothing
        end
        i = i+1
    end
    println("$(round(sum(AA_vector) - ((n_end - n_start - 1) * WATER_MW); digits = 3)) g/mol")
    println("$(round(((sum(AA_vector) - ((n_end - n_start - 1) * WATER_MW))/1000); digits = 3)) kDa")
end
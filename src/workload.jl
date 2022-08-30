using BioSequences
using SnoopPrecompile

# BioSequences define a whole bunch of types and methods, most of which is never used in any workload.
# This workload here uses what I consider to be the most common operations,
# on the most common types.
# This is intended to strike a balance between precompilign that which the user probably needs, while
# not wasting time loading code the user will never use.
# The code not cached here will have to be precompiled in downstream packages

@precompile_all_calls begin
    seqs = [
        aa"TAGCW"
        dna"ATCGA"
    ]
    for seq in [seqs; map(i -> view(i, 1:5), seqs)]
        # printing
        String(seq)
        print(IOBuffer(), seq)

        hash(seq)

        # indexing
        seq[1]
        seq[2:3]
        seq[2] = seq[3]
        seq[2:3] = seq[3:4]

        # join
        join([seq, seq])
        join((first(seq), last(seq)))

        # pred
        hasambiguity(seq)

        # transformations
        reverse(seq)
        ungap(seq)

        q1 = ExactSearchQuery(seq[1:2])
        findfirst(q1, seq)
        findlast(q1, seq)
        findall(q1, seq)
        occursin(q1, seq)

        q2 = ApproximateSearchQuery(seq[2:4])
        findfirst(q2, 1, seq)
        findlast(q2, 1, seq)
        occursin(q2, 1, seq)
    end

    for seq in seqs
        seq[[true, false, true, true, false]]
        seq[collect(2:3)]
        ungap!(seq)
    end

    # Nucleotide
    for seq in seqs[2:2]
        ispalindromic(seq)
        iscanonical(seq)

        canonical(seq)
        reverse_complement!(seq)
        reverse_complement(seq)
        complement(seq)
        translate(seq[1:3])
        gc_content(seq)
    end

    # Random
    randdnaseq(5)
    randrnaseq(5)
    randaaseq(5)
    randseq(DNAAlphabet{2}(), 10)    
end

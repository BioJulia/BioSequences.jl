###
### Bitparallel 
###

"""
A generator for efficient bitwise sequence operations.
"""
function compile_bitpar(funcname::Symbol;
                        arguments::Tuple  =  (),
                        init_code::Expr   = :(),
                        head_code::Expr   = :(),
                        body_code::Expr   = :(),
                        tail_code::Expr   = :(),
                        return_code::Expr = :())
    functioncode = :(function $(funcname)() end)
    for arg in arguments
        push!(functioncode.args[1].args, arg)
    end
    functioncode.args[2] = quote
        $(init_code)
        ind = bitindex(seq, 1)
        stop = bitindex(seq, lastindex(seq) + 1)
        data = seq.data
        @inbounds begin
            if !iszero(offset(ind)) & (ind < stop)
                # align the bit index to the beginning of a block boundary
                o = offset(ind)
                mask = bitmask(stop - ind)
                n_bits_masked = ifelse(index(stop) == index(ind), count_zeros(mask), o)
                chunk = (data[index(ind)] >> o) & mask
                $(head_code)
                ind += 64 - o
            end
            
            lastind = index(stop - bits_per_symbol(Alphabet(seq)))
            lastind -= !iszero(offset(stop))
            for i in index(ind):lastind
                chunk = data[i]
                $(body_code)
                ind += 64
            end
            
            if ind < stop
                n_bits_masked = 64 - offset(stop)
                chunk = data[index(ind)] & bitmask(64 - n_bits_masked)
                $(tail_code)
            end
        end
        $(return_code)
    end
    return functioncode
end

function get_arg_name(arg)
    if typeof(arg) === Expr
        if arg.head === :(::)
            return first(arg.args)
        end
    end
    return arg
end

function check_arguments(args::Tuple)
    # Check all arguments are symbols or expressions
    for arg in args
        if !(isa(arg, Symbol) || isa(arg, Expr))
            error("Argument ", arg, " is not an expression or symbol")
        end
    end
    
    a1 = first(args)
    if isa(a1, Symbol)
        @assert a1 == :seqa
    elseif isa(a1, Expr)
        @assert a1.head === :(::)
        @assert first(a1.args) === :seqa
    end
    
    a2 = args[2]
    if isa(a2, Symbol)
        @assert a2 == :seqb
    elseif isa(a2, Expr)
        @assert a2.head === :(::)
        @assert first(a2.args) === :seqb
    end
end

function add_arguments!(exp, args)
    check_arguments(args)
    @assert exp.head === :function
    if exp.args[1].head === :call
        for arg in args
            push!(exp.args[1].args, arg)
        end
    elseif exp.args[1].head === :where
        @assert exp.args[1].args[1].head === :call
        for arg in args
            push!(exp.args[1].args[1].args, arg)
        end
    else
        error("Expression in function is not a :call or :where!")
    end
end

function compile_2seq_bitpar(funcname::Symbol;
                             arguments::Tuple  =  (),
                             parameters::Tuple =  (),
                             init_code::Expr   = :(),
                             head_code::Expr   = :(),
                             body_code::Expr   = :(),
                             tail_code::Expr   = :(),
                             return_code::Expr = :())
    
    # TODO: Check the parameters provided.
    
    if isempty(parameters)
        functioncode = :(function $(funcname)() end)
    else
        functioncode = :(function $(funcname)() where {$(parameters...)} end)
    end
    
    add_arguments!(functioncode, arguments)
    
    argument_names = map(get_arg_name, arguments)
    argument_names = tuple(argument_names[2], argument_names[1], argument_names[3:end]...)
    
    functioncode.args[2] = quote
        if length(seqa) > length(seqb)
            return $funcname($(argument_names...))
        end
        @assert length(seqa) ≤ length(seqb)
        
        nexta = bitindex(seqa, 1)
        stopa = bitindex(seqa, lastindex(seqa) + 1)
        nextb = bitindex(seqb, 1)
        stopb = bitindex(seqb, lastindex(seqb) + 1)
        adata = seqa.data
        bdata = seqb.data
        
        $(init_code)
        
        # The first thing we need to sort out is to correctly align the head of
        # sequence / subsequence `a`s data is aligned such that the offset of
        # `nexta` is essentially reduced to 0.
        # With sequence / subsequence `a` aligned, from there, we only need to
        # worry about the alignment of sequence / subsequence `b` with respect
        # to `a`.
        if nexta < stopa && offset(nexta) != 0
            # Here we shift the first data chunks to the right so as the first
            # nucleotide of the seq/subseq is the first nibble / pair of bits.
            x = adata[index(nexta)] >> offset(nexta)
            y = bdata[index(nextb)] >> offset(nextb)
            # Here it was assumed that there is something to go and get from
            # the next chunk of `b`, yet that may not be true.
            # We know that if this is not true of `b`, then it is certainly not
            # true of `a`.
            # We check if the end of the sequence is not contained in the same
            # integer like so: `64 - offset(nextb) < stopb - nextb`.
            #
            # This edge case was found and accounted for by Ben J. Ward @BenJWard.
            # Ask this maintainer for more information.
            if offset(nextb) > offset(nexta) && 64 - offset(nextb) < stopb - nextb
                y |= bdata[index(nextb) + 1] << (64 - offset(nextb))
            end
            # Here we need to check something, we need to check if the
            # integer of `a` we are currently aligning contains the end of
            # seq/subseq `a`. Because if so it's something we need to take into
            # account of when we mask x and y.
            #
            # In other words if `64 - offset(nexta) > stopa - nexta, we know
            # seq or subseq a's data ends before the end of this data chunk,
            # and so the mask used needs to be defined to account for this:
            # `mask(stopa - nexta)`, otherwise the mask simply needs to be
            # `mask(64 - offset(nexta))`.
            #
            # This edge case was found and accounted for by Ben Ward @Ward9250.
            # Ask this maintainer for more information.
            offs = ifelse(64 - offset(nexta) > stopa - nexta, stopa - nexta, 64 - offset(nexta))
            m = bitmask(offs)
            
            x &= m
            y &= m
            
            $(head_code)
            
            # Here we move our current position markers by k, meaning they move
            # to either, A). The next integer, or B). The end of the sequence if
            # it is in the current integer.
            nexta += offs
            nextb += offs
        end
        
        if offset(nextb) == 0  # data are aligned with each other
            while stopa - nexta ≥ 64 # Iterate through body of data
                x = adata[index(nexta)]
                y = bdata[index(nextb)]
                
                $(body_code)
                
                nexta += 64
                nextb += 64
            end

            if nexta < stopa
                x = adata[index(nexta)]
                y = bdata[index(nextb)]

                offs = stopa - nexta
                m = bitmask(offs)
                x &= m
                y &= m
                
                $(tail_code)
            end
        elseif nexta < stopa # Data are unaligned
            y = bdata[index(nextb)]
            nextb += 64
            # Note that here, updating `nextb` by 64, increases the chunk index,
            # but the `offset(nextb)` will remain the same.
            while stopa - nexta ≥ 64 # processing body of data
                x = adata[index(nexta)]
                z = bdata[index(nextb)]
                y = y >> offset(nextb) | z << (64 - offset(nextb))
                
                $(body_code)

                y = z
                nexta += 64
                nextb += 64
            end
            
            if nexta < stopa # processing tail of data
                x = adata[index(nexta)]
                y = y >> offset(nextb)
                if 64 - offset(nextb) < stopa - nexta
                    y |= bdata[index(nextb)] << (64 - offset(nextb))
                end
                offs = stopa - nexta
                m = bitmask(offs)
                x &= m
                y &= m
                
                $(tail_code)
            end
        end
        $(return_code)
    end
    return functioncode
end

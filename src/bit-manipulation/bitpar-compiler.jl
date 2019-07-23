###
### Bitparallel 
###

"""
A generator for efficient bitwise sequence operations.
"""
function compile_bitpar(funcname::Symbol;
                        arguments::Tuple  = (),
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
        i = bitindex(seq, 1)
        stop = bitindex(seq, lastindex(seq) + 1)
        @inbounds begin
            if offset(i) != 0 && i < stop
                # align the bit index to the beginning of a block boundary
                o = offset(i)
                chunk = (seq.data[index(i)] >> o) & bitmask(stop - i)
                $(head_code)
                i += 64 - o
                @assert offset(i) == 0
            end
            while i ≤ stop - 64
                chunk = seq.data[index(i)]
                $(body_code)
                i += 64
            end
            if i < stop
                chunk = seq.data[index(i)] & bitmask(offset(stop))
                $(tail_code)
            end
        end
        $(return_code)
    end
    return functioncode
end
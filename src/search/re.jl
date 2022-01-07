# Regular Expression
# ==================
#
# Regular expression sequence search tools.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# String Decorators
# -----------------

"""
    biore"PAT"sym

Construct a PCRE `BioRegex` from pattern `PAT`. `sym` can be any of `d`/`dna`, `r`/`rna`,
or `a`/`aa`. BioRegex can also be constructed directly using e.g. `BioRegex{DNA}("[TA]G")`.
"""
macro biore_str(pat, opt...)
    if isempty(opt)
        error("symbol option is required: d(na), r(na), or a(a)")
    end
    @assert length(opt) == 1
    opt = opt[1]

    if opt ∈ ("d", "dna")
        :($(RE.Regex{DNA}(pat, :pcre)))
    elseif opt ∈ ("r", "rna")
        :($(RE.Regex{RNA}(pat, :pcre)))
    elseif opt ∈ ("a", "aa")
        :($(RE.Regex{AminoAcid}(pat, :pcre)))
    else
        error("invalid symbol option: $(opt)")
    end
end

"""
    prosite"PAT"

Equivalent to `BioRegex{AminoAcid}("PAT", syntax=:prosite)`, but evaluated at parse-time.
"""
macro prosite_str(pat)
    :($(RE.Regex{AminoAcid}(pat, :prosite)))
end


module RE

import BioSequences

# Syntax tree
# -----------

mutable struct SyntaxTree
    head::Symbol
    args::Vector{Any}

    function SyntaxTree(head, args)
        @assert head ∈ (
            :|,
            :*, :+, Symbol("?"), :range,
            Symbol("*?"), Symbol("+?"), Symbol("??"), Symbol("range?"),
            :set, :compset, :sym, :bits,
            :capture, :nocapture, :concat, :head, :last)
        return new(head, args)
    end
end

expr(head, args) = SyntaxTree(head, args)

function charset(s)
    set = Set{Char}()
    union!(set, s)
    union!(set, lowercase(s))
    return set
end

# list of symbols available for each symbol type
const symbols = IdDict{Any, Set{Char}}(
    BioSequences.DNA => charset("ACGTMRWSYKVHDBN"),
    BioSequences.RNA => charset("ACGUMRWSYKVHDBN"),
    BioSequences.AminoAcid     => charset("ARNDCQEGHILKMFPSTWYVOUBJZX"))

macro check(ex, err)
    esc(quote
        if !$ex
            throw($err)
        end
    end)
end

function parse(::Type{T}, pat::AbstractString) where {T}
    parens = Char[]  # stack of parens
    ex, _ = parserec(T, pat, iterate(pat), parens)
    @check isempty(parens) ArgumentError("'(' is not closed")
    return ex
end

function parserec(::Type{T}, pat, s, parens) where {T}
    args = []
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if c == '*'
            @check !isempty(args) ArgumentError("unexpected '*'")
            arg = pop!(args)
            push!(args, expr(:*, [arg]))
        elseif c == '+'
            @check !isempty(args) ArgumentError("unexpected '+'")
            arg = pop!(args)
            push!(args, expr(:+, [arg]))
        elseif c == '?'
            @check !isempty(args) ArgumentError("unexpected '?'")
            arg = pop!(args)
            if arg.head ∈ (:*, :+, Symbol("?"), :range)
                # lazy quantifier
                push!(args, expr(Symbol(arg.head, '?'), arg.args))
            else
                # zero-or-one quantifier
                push!(args, expr(Symbol("?"), [arg]))
            end
        elseif c == '{'
            @check !isempty(args) ArgumentError("unexpected '{'")
            rng, s = parserange(pat , s)
            arg = pop!(args)
            push!(args, expr(:range, [rng, arg]))
        elseif c == '|'
            arg1 = expr(:concat, args)
            arg2, s = parserec(T, pat, s, parens)
            args = []
            push!(args, expr(:|, [arg1, arg2]))
        elseif c == '['
            setexpr, s = parseset(T, pat, s)
            push!(args, setexpr)
        elseif c == '('
            push!(parens, '(')
            head = :capture
            if peek(pat, s) == '?'
                _, state = s
                s = iterate(pat, state)
                if peek(pat, s) == ':'
                    head = :nocapture
                    _, state = s
                    s = iterate(pat, state)
                end
            end
            arg, s = parserec(T, pat, s, parens)
            push!(args, expr(head, [arg]))
        elseif c == ')'
            @check !isempty(parens) && parens[end] == '(' ArgumentError("unexpected ')'")
            pop!(parens)
            return expr(:concat, args), s
        elseif c == '^'
            push!(args, expr(:head, []))
        elseif c == '$'
            push!(args, expr(:last, []))
        elseif isspace(c)
            # skip
        elseif c ∈ symbols[T]
            push!(args, expr(:sym, [convert(T, c)]))
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(:concat, args), s
end

function parserange(pat, s)
    lo = hi = -1
    comma = false
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if isdigit(c)
            d = c - '0'
            if comma
                if hi < 0
                    hi = 0
                end
                hi = 10hi + d
            else
                if lo < 0
                    lo = 0
                end
                lo = 10lo + d
            end
        elseif c == ','
            comma = true
        elseif c == '}'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    if comma
        if lo < 0 && hi < 0
            throw(ArgumentError("invalid range"))
        elseif lo ≥ 0 && hi < 0
            return (lo,), s
        elseif lo < 0 && hi ≥ 0
            return (0, hi), s
        else  # lo ≥ 0 && hi ≥ 0
            if lo > hi
                throw(ArgumentError("invalid range"))
            end
            return (lo, hi), s
        end
    else
        return lo, s
    end
end

function peek(pat, s)
    if s === nothing
        throw(ArgumentError("unexpected end of pattern"))
    end
    return s[1]
end

function parseset(::Type{T}, pat, s) where {T}
    if peek(pat, s) == '^'
        head = :compset
        _, state = s
        s = iterate(pat, state)
    else
        head = :set
    end
    set = T[]
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if c ∈ symbols[T]
            push!(set, convert(T, c))
        elseif c == ']'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(head, set), s
end

function parse_prosite(pat)
    s = iterate(pat)
    args = []
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if c == '['
            set, s = parseset_prosite(pat, s, ']')
            push!(args, set)
        elseif c == '{'
            set, s = parseset_prosite(pat, s, '}')
            push!(args, set)
        elseif c == '('
            @check !isempty(args) ArgumentError("unexpected '('")
            rng, s = parserange_prosite(pat, s)
            arg = pop!(args)
            push!(args, expr(:range, [rng, arg]))
        elseif c == '-'
            # concat
            continue
        elseif c == 'x'
            push!(args, expr(:sym, [BioSequences.AA_X]))
        elseif c == '<'
            push!(args, expr(:head, []))
        elseif c == '>'
            push!(args, expr(:last, []))
        elseif c ∈ symbols[BioSequences.AminoAcid]
            push!(args, expr(:sym, [convert(BioSequences.AminoAcid, c)]))
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    return expr(:concat, args)
end

function parserange_prosite(pat, s)
    lo = hi = -1
    comma = false
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if isdigit(c)
            d = c - '0'
            if comma
                if hi < 0
                    hi = 0
                end
                hi = 10hi + d
            else
                if lo < 0
                    lo = 0
                end
                lo = 10lo + d
            end
        elseif c == ','
            comma = true
        elseif c == ')'
            break
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
    if comma
        if lo < 0 || hi < 0
            throw(ArgumentError("invalid range"))
        end
        return (lo, hi), s
    else
        if lo < 0
            throw(ArgumentError("invalid range"))
        end
        return lo, s
    end
end

function parseset_prosite(pat, s, close)
    set = BioSequences.AminoAcid[]
    while s !== nothing
        c, state = s
        s = iterate(pat, state)
        if c ∈ symbols[BioSequences.AminoAcid]
            push!(set, convert(BioSequences.AminoAcid, c))
        elseif c == close
            if close == ']'
                return expr(:set, set), s
            elseif close == '}'
                return expr(:compset, set), s
            end
            @assert false
        else
            throw(ArgumentError("unexpected input: '$c'"))
        end
    end
end

function bits2sym(::Type{T}, bits::UInt32) where {T}
    for x in BioSequences.alphabet(T)
        if BioSequences.compatbits(x) == bits
            return x
        end
    end
    error("bits are not found")
end

mask(::Type{T}) where {T<:BioSequences.NucleicAcid} = (UInt32(1) << 4) - one(UInt32)
@assert Int(reinterpret(UInt8, BioSequences.AA_U)) + 1 == 22  # check there are 22 unambiguous amino acids
mask(::Type{BioSequences.AminoAcid}) = (UInt32(1) << 22) - one(UInt32)

function desugar(::Type{T}, tree::SyntaxTree) where {T}
    head = tree.head
    args = tree.args
    if head == :+
        # e+ => ee*
        head = :concat
        args = [args[1], expr(:*, [args[1]])]
    elseif head == Symbol("+?")
        # e+? => ee*?
        head = :concat
        args = [args[1], expr(Symbol("*?"), [args[1]])]
    elseif head == Symbol("?")
        # e? => e|
        head = :|
        args = [args[1], expr(:concat, [])]
    elseif head == Symbol("??")
        # e?? => |e
        head = :|
        args = [expr(:concat, []), args[1]]
    elseif head == :sym
        head = :bits
        args = [BioSequences.compatbits(args[1])]
    elseif head == :set
        bits = UInt32(0)
        for arg in args
            bits |= BioSequences.compatbits(arg)
        end
        head = :bits
        args = [bits]
    elseif head == :compset
        bits = UInt32(0)
        for arg in args
            bits |= BioSequences.compatbits(arg)
        end
        head = :bits
        args = [~bits & mask(T)]
    elseif head == :nocapture
        head = :concat
    elseif head == :range || head == Symbol("range?")
        rng = args[1]
        pat = args[2]
        greedy = head == :range
        if isa(rng, Int)
            # e{m} => eee...e
            #         |<-m->|
            head = :concat
            args = fill(pat, rng)
        elseif isa(rng, Tuple{Int})
            # e{m,} => eee...ee*   (greedy)
            #          |<-m->|
            # e{m,} => eee...ee*?  (lazy)
            #          |<-m->|
            head = :concat
            args = fill(pat, rng[1])
            if greedy
                push!(args, expr(:*, [pat]))
            else
                push!(args, expr(Symbol("*?"), [pat]))
            end
        elseif isa(rng, Tuple{Int,Int})
            # e{m,n} => eee...eee|eee...ee|...|eee...e  (greedy)
            #           |<- n ->|              |<-m->|
            # e{m,n} => eee...e|eee...ee|...|eee...eee  (lazy)
            #           |<-m->|              |<- n ->|
            m, n = rng
            @assert m ≤ n
            if m == n
                head = :concat
                args = fill(pat, m)
            else
                head = :|
                f = (k, t) -> expr(:|, [expr(:concat, fill(pat, k)), t])
                if greedy
                    args = foldr(f, n:-1:m+1, init=expr(:concat, fill(pat, m))).args
                else
                    args = foldr(f, m:1:n-1, init=expr(:concat, fill(pat, n))).args
                end
            end
        else
            @assert false "invalid AST"
        end
    end
    i = 1
    for i in 1:lastindex(args)
        args[i] = desugar(T, args[i])
    end
    return expr(head, args)
end

function desugar(::Type{T}, atom) where {T}
    return atom
end


# Compiler
# --------

#  tag  | operation |           meaning
# ------|-----------|----------------------------
# 0b000 | match     | the pattern matches
# 0b001 | bits b    | matches b in bit-wise way
# 0b010 | jump l    | jump to l
# 0b011 | push l    | push l and go to next
# 0b100 | save i    | save string state to i
# 0b101 | head      | matches the head of string
# 0b110 | last      | matches the last of string
# 0b111 | fork l    | push next and go to l

primitive type Op 32 end

const MatchTag = UInt32(0b000) << 29
const BitsTag  = UInt32(0b001) << 29
const JumpTag  = UInt32(0b010) << 29
const PushTag  = UInt32(0b011) << 29
const SaveTag  = UInt32(0b100) << 29
const HeadTag  = UInt32(0b101) << 29
const LastTag  = UInt32(0b110) << 29
const ForkTag  = UInt32(0b111) << 29

# constructors
match() = reinterpret(Op, MatchTag)
bits(b::UInt32) = reinterpret(Op, BitsTag | b)
jump(l::Int) = reinterpret(Op, JumpTag | UInt32(l))
push(l::Int) = reinterpret(Op, PushTag | UInt32(l))
save(l::Int) = reinterpret(Op, SaveTag | UInt32(l))
head() = reinterpret(Op, HeadTag)
last() = reinterpret(Op, LastTag)
fork(l::Int) = reinterpret(Op, ForkTag | UInt32(l))

const operand_mask = (UInt32(1) << 29) - one(UInt32)

function tag(op::Op)
    return reinterpret(UInt32, op) & ~operand_mask
end

function operand(op::Op)
    return reinterpret(UInt32, op) &  operand_mask
end

function Base.show(io::IO, op::Op)
    t = tag(op)
    x = operand(op)
    if t == MatchTag
        print(io, "match")
    elseif t == BitsTag
        print(io, "bits ", repr(x))
    elseif t == JumpTag
        print(io, "jump ", x)
    elseif t == PushTag
        print(io, "push ", x)
    elseif t == SaveTag
        print(io, "save ", x)
    elseif t == HeadTag
        print(io, "head")
    elseif t == LastTag
        print(io, "last")
    elseif t == ForkTag
        print(io, "fork ", x)
    else
        @assert false
    end
end

function print_program(prog::Vector{Op})
    L = lastindex(prog)
    for l in 1:L
        print(lpad(l, ndigits(L)), ": ", prog[l])
        if l != lastindex(prog)
            println()
        end
    end
end

function compile(tree::SyntaxTree)
    code = Op[]
    push!(code, save(1))
    compilerec!(code, tree, 2)
    push!(code, save(2))
    push!(code, match())
    return code
end

function compilerec!(code, tree::SyntaxTree, k)
    h = tree.head
    args = tree.args
    if h == :bits
        push!(code, bits(UInt32(args[1])))
    elseif h == :concat
        for arg in args
            k = compilerec!(code, arg, k)
        end
    elseif h == :*
        push!(code, push(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(l))
        code[l] = push(length(code) + 1)
    elseif h == Symbol("*?")
        push!(code, fork(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(l))
        code[l] = fork(length(code) + 1)
    elseif h == :|
        push!(code, push(0))  # placeholder
        l = length(code)
        k = compilerec!(code, args[1], k)
        push!(code, jump(0))  # placeholder
        code[l] = push(length(code) + 1)
        l = length(code)
        k = compilerec!(code, args[2], k)
        code[l] = jump(length(code) + 1)
    elseif h == :capture
        k′ = k
        push!(code, save(2k′-1))
        k = compilerec!(code, args[1], k + 1)
        push!(code, save(2k′))
    elseif h == :head
        push!(code, head())
    elseif h == :last
        push!(code, last())
    else
        @assert false "invalid tree"
    end
    return k
end


# Virtual machine
# ---------------

"""
    Regex{T}(pattern::AbstractString, syntax=:pcre)

Biological regular expression to seatch for `pattern` in sequences of type `T`, where
`T` can be `DNA`, `RNA`, and `AminoAcid`. `syntax` can be `:pcre` or `:prosite` for AminoAcid
acids.
"""
struct Regex{T}
    pat::String       # regular expression pattern (for printing)
    code::Vector{Op}  # compiled code
    nsaves::Int       # the number of `save` operations in `code`

    function Regex{T}(pat::AbstractString, syntax=:pcre) where T
        if syntax == :pcre
            ast = desugar(T, parse(T, pat))
        elseif syntax == :prosite
            if T != BioSequences.AminoAcid
                throw(ArgumentError("alphabet must be AminoAcid for PROSITE syntax"))
            end
            ast = desugar(BioSequences.AminoAcid, parse_prosite(pat))
        else
            throw(ArgumentError("invalid syntax: $syntax"))
        end
        code = compile(ast)
        nsaves = 0
        for op in code
            nsaves += tag(op) == SaveTag
        end
        @assert iseven(nsaves)
        return new{T}(pat, code, nsaves)
    end
end

function Base.show(io::IO, re::Regex{T}) where {T}
    if T == BioSequences.DNA
        opt = "dna"
    elseif T == BioSequences.RNA
        opt = "rna"
    elseif T == BioSequences.AminoAcid
        opt = "aa"
    else
        assert(false)
    end
    print(io, "biore\"", re.pat, "\"", opt)
end

# This is useful when testing
function Base.:(==)(x::Regex, y::Regex)
    typeof(x) == typeof(y) && x.pat == y.pat && x.code == y.code && x.nsaves == y.nsaves
end

"""
    RegexMatch

Result of matching by `Regex`.

julia> match(biore"A(C[TG])+N(CA)"d, dna"ACGACA")
RegexMatch("ACGACA", 1="CG", 2="", 3="CA")
"""
struct RegexMatch{S}
    seq::S
    captured::Vector{Int}
end

function Base.show(io::IO, m::RegexMatch)
    print(io, "RegexMatch(")
    for k in 1:div(length(m.captured), 2)
        if k > 1
            print(io, ", ", k - 1, '=')
        end
        print(io, '"', m.seq[m.captured[2k-1]:m.captured[2k]-1], '"')
    end
    print(io, ')')
end

"""
    matched(match::BioRegexMatch)

Return the matched pattern as a `BioSequence`.
"""
function matched(m::RegexMatch{S}) where {S}
    return m.seq[m.captured[1]:m.captured[2]-1]
end


"""
    captured(match::BioRegexMatch)

Return a vector of the captured patterns, where a pattern not captured is `nothing`.

# Examples:
```
julia> captured(biore"A(C[TG])+N"d, dna"ACGAA"")
2-element Vector{Union{Nothing, LongDNA{4}, LongNuc{4, DNAAlphabet{4}}}}:
 CG
 nothing
```
"""
function captured(m::RegexMatch{S}) where {S}
    return [m.captured[2k-1] != 0 && m.captured[2k] != 0 ?
            m.seq[m.captured[2k-1]:m.captured[2k]-1] : nothing
            for k in 2:div(length(m.captured), 2)]
end

function checkeltype(re::Regex{T}, seq::BioSequences.BioSequence) where {T}
    if eltype(seq) != T
        throw(ArgumentError("element type of sequence doesn't match with regex"))
    end
end

function Base.match(re::Regex{T}, seq::BioSequences.BioSequence, start::Integer=1) where {T}
    checkeltype(re, seq)

    # use the first unambiguous symbol in the regular expression to find the
    # next starting position of pattern matching; this improves the performance
    if length(re.code) ≥ 2 && tag(re.code[2]) == BitsTag && count_ones(operand(re.code[2])) == 1
        firstsym = bits2sym(T, operand(re.code[2]))
    else
        firstsym = BioSequences.gap(T)
    end

    # a thread is `(<program counter>, <sequence's iterator state>)`
    threads = Stack{Tuple{Int,Int}}()
    captured = Vector{Int}(undef, re.nsaves)
    s = start
    while true
        if firstsym != BioSequences.gap(T)
            s = findnext(isequal(firstsym), seq, s)
            if s === nothing
                break
            end
        end
        empty!(threads)
        push!(threads, (1, s))
        fill!(captured, 0)
        if runmatch!(threads, captured, re, seq)
            return RegexMatch(seq, captured)
        end
        s += 1
        if s > lastindex(seq)
            break
        end
    end
    return nothing
end

function Base.findfirst(re::Regex{T}, seq::BioSequences.BioSequence, start::Integer=1) where {T}
    checkeltype(re, seq)
    m = Base.match(re, seq, start)
    if m === nothing
        return nothing
    else
        return m.captured[1]:m.captured[2]-1
    end
end

struct RegexMatchIterator{T,S}
    re::Regex{T}
    seq::S
    overlap::Bool

    function RegexMatchIterator{T,S}(re::Regex{T}, seq::S, overlap::Bool) where {T,S}
        checkeltype(re, seq)
        return new{T,S}(re, seq, overlap)
    end
end

function Base.IteratorSize(::Type{RegexMatchIterator{T,S}}) where {T,S}
    return Base.SizeUnknown()
end

function Base.eltype(::Type{RegexMatchIterator{T,S}}) where {T,S}
    return RegexMatch{S}
end

function Base.iterate(iter::RegexMatchIterator)
    threads = Stack{Tuple{Int,Int}}()
    captured = Vector{Int}(undef, iter.re.nsaves)
    s = 1
    push!(threads, (1, s))
    fill!(captured, 0)
    state = advance!(threads, captured, s, iter.re, iter.seq, iter.overlap)
    return iterate(iter, state)
end

function Base.iterate(iter::RegexMatchIterator, state)
    item, threads, captured, s = state
    if item === nothing
        return nothing
    else
        item, advance!(threads, captured, s, iter.re, iter.seq, iter.overlap)
    end
end

function advance!(threads, captured, s, re, seq, overlap)
    while true
        if runmatch!(threads, captured, re, seq)
            if !overlap
                empty!(threads)
                s = captured[2]
                if s <= lastindex(seq)
                    push!(threads, (1, s))
                end
            end
            return RegexMatch(seq, copy(captured)), threads, captured, s
        end
        if !isempty(seq)
            s += 1
        end
        if s > lastindex(seq)
            break
        end
        push!(threads, (1, s))
        fill!(captured, 0)
    end
    return nothing, threads, captured, s
end

function Base.eachmatch(re::Regex{T}, seq::BioSequences.BioSequence, overlap::Bool = true) where {T}
    checkeltype(re, seq)
    return RegexMatchIterator{T,typeof(seq)}(re, seq, overlap)
end

function Base.occursin(re::Regex{T}, seq::BioSequences.BioSequence) where {T}
    return Base.match(re, seq) !== nothing
end

# simple stack
mutable struct Stack{T}
    top::Int
    data::Vector{T}

    function Stack{T}(sz::Int=0) where T
        return new{T}(0, Vector{T}(undef, sz))
    end
end

function Base.isempty(stack::Stack)
    return stack.top == 0
end

@inline function Base.push!(stack::Stack{T}, x::T) where {T}
    if stack.top + 1 > length(stack.data)
        push!(stack.data, x)
    else
        stack.data[stack.top+1] = x
    end
    stack.top += 1
    return stack
end

@inline function Base.pop!(stack::Stack)
    item = stack.data[stack.top]
    stack.top -= 1
    return item
end

@inline function Base.empty!(stack::Stack)
    stack.top = 0
    return stack
end

# run pattern maching using the current `threads` stack. Captured positions are
# stored in the preallocated `captured`. If match is found, this returns `true`
# immediately; otherwise returns `false`.
function runmatch!(threads::Stack{Tuple{Int,Int}},
                   captured::Vector{Int},
                   re::Regex, seq::BioSequences.BioSequence)
    while !isempty(threads)
        pc::Int, s = pop!(threads)
        while true
            op = re.code[pc]
            t = tag(op)
            if t == BitsTag
                if s > lastindex(seq)
                    break
                end
                sym = seq[s]
                s += 1
                if BioSequences.compatbits(sym) & operand(op) != 0
                    pc += 1
                else
                    break
                end
            elseif t == JumpTag
                pc = operand(op)
            elseif t == PushTag
                push!(threads, (convert(Int, operand(op)), s))
                pc += 1
            elseif t == ForkTag
                push!(threads, (pc + 1, s))
                pc = operand(op)
            elseif t == SaveTag
                captured[operand(op)] = s
                pc += 1
            elseif t == HeadTag
                if s == 1
                    pc += 1
                else
                    break
                end
            elseif t == LastTag
                if s == lastindex(seq) + 1
                    pc += 1
                else
                    break
                end
            elseif t == MatchTag
                return true
            end
        end
    end
    return false
end

end  # module RE

# exported from BioSequences
const matched = RE.matched
const captured = RE.captured
const BioRegex = RE.Regex
const BioRegexMatch = RE.RegexMatch

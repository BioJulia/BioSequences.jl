# * iscertain + !iscertain
#   Four bit: Certain bitcount OR enumerate_nibbles & 0x1111111 + tz
#   AA: 0x00-0x15 or 0x1a
#

# * isambiguous + !isambiguous
#   Four bit: Ambiguous bitcount OR enumerate_nibbles & 0xEEE.. + tz
#   AA: 0x16-0x19

# Make 0 for any other bitpattern than 0x1, 0x2, 0x4, 0x8.

# Zeros out all nibbles except the ones with 4bit A,C,G or T
function iscertain_kernel(x::UInt64)
    x = enumerate_nibbles(x)
    y =  (x & 0x4444444444444444) >> 2
    y |= (x & 0x2222222222222222) >> 1
    x & ~y & 0x1111111111111111
end

# Zeros out all which is not a certain AA
# TODO: This is not efficient.
function iscertain_aa_kernel(x::UInt64)
    y = reinterpret(NTuple{8, UInt8}, x)
    y = map(i -> ((i < 0x16) | (i == 0x1a)), y)
    reinterpret(UInt64, y)
end


# Zeros out all nibbles encoding 4bit A,C,G or T
uncertain_kernel(x::UInt64) = enumerate_nibbles(x) ⊻ 0x1111111111111111

# Zeros out A, C, G, T or Gap
function ambiguous_kernel(x::UInt64)
    # The y part makes every nibble 0xF, unless it's 0 to begin with
    y = x | (x >>> 2)
    y |= y >>> 1
    y &= 0x1111111111111111
    y *= 0xF
    uncertain_kernel(x) & y
end

# Zeros out all except A, C, G, T and Gap
function unambiguous_kernel(x::UInt64)
    y = enumerate_nibbles(x)
    y = (y >>> 1) | (y >>> 2)
    y &= 0x1111111111111111
    y ⊻= 0x1111111111111111
end

# For debug
#=
function make_all_nibbles()
    x = UInt64(0)
    for i in 0:15
        x = (x << 4) | i
    end
    x
end
=#
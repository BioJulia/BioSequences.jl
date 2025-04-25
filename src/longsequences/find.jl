# * iscertain + !iscertain
#   Four bit: Certain bitcount OR enumerate_nibbles & 0x1111111 + tz
#   AA: 0x00-0x15 or 0x1a
#

# * isambiguous + !isambiguous
#   Four bit: Ambiguous bitcount OR enumerate_nibbles & 0xEEE.. + tz
#   AA: 0x16-0x19

# * isgap + !isgap
#   Four: tz
#   AA: 0x1b

# * isGC + !isGC
# * isequal / == 

function findnext(iscertain, seq::SeqOrView{<:NucleicAcidAlphabet{4}}, from::Integer)
    from = Int(from)::Int
    @boundscheck (from < 1 && throw(BoundsError(seq, from)))
end
```@meta
CurrentModule = BioSequences
DocTestSetup = quote
    using BioSequences
end
```
# BitIndex

The BitIndex type is an important type in BioSequences internally.
BioSequence types store their elements in a vector in a succinct and packed
form, allowing multiple symbols (dna, rna, amino-acids, etc.) to be stored
inside a single word. This has great performance benefits but makes some things
more tricky. For example iterating though the encoded data of the sequence or
keeping track of where in the binary array the elements begin and end.

We simplify this by keeping the length of encoded binary bits in a sequence fixed.
Hence a character at arbitrary position can be extracted in a constant time.
We then use this BitIndex type which represents the position of an element in the
encoded bits as an index, and an offset. I.e Which chunk an element is stored in,
and its offset inside of that chunk.

You can think of this visually:

#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

In the above diagram, the element pointed to is marked by X, `i` is the BitIndex
and in this case the chunk size is 64 bits (although the BitIndex type supports
other chunk sizes).

The BitIndex type has two parameters: `N` which is the number of bits used for
each element, and `W`, which is the word (or chunk) size.
 
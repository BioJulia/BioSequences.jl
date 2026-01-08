@testset "Bit manipulation" begin
    function test_bit_reverse(u, bps)
        rv = BioSequences.reversebits(u, BioSequences.BitsPerSymbol{bps}())
        ustr = string(u; base=2, pad=8 * sizeof(u))
        rvs = join(reverse!(collect(Iterators.partition(ustr, bps))))
        @test rv == parse(typeof(u), rvs; base=2)
    end

    # Test UInt128 with various BitsPerSymbol values (powers of two: 1, 2, 4, 8, 16, 32, 64)
    for bps in [1, 2, 4, 8, 16, 32, 64]
        # Edge cases
        test_bit_reverse(UInt128(0), bps)
        test_bit_reverse(typemax(UInt128), bps)

        # Small values
        test_bit_reverse(UInt128(1), bps)
        test_bit_reverse(UInt128(0xff), bps)
        test_bit_reverse(UInt128(0xffff), bps)

        # # Medium values
        test_bit_reverse(UInt128(0x123456789abcdef0), bps)
        test_bit_reverse(UInt128(0xfedcba9876543210), bps)

        # # Large values
        test_bit_reverse(UInt128(0xb1d318f6d8b882f1ee180f1bcd8f8727), bps)
        test_bit_reverse(UInt128(0xffffffffffffffff0000000000000000), bps)
        test_bit_reverse(UInt128(0x0000000000000000ffffffffffffffff), bps)
        test_bit_reverse(UInt128(0xaaaaaaaaaaaaaaaa5555555555555555), bps)
        test_bit_reverse(UInt128(0x5555555555555555aaaaaaaaaaaaaaaa), bps)

        # # Random-looking patterns
        test_bit_reverse(UInt128(0xdeadbeefcafebabe0123456789abcdef), bps)
        test_bit_reverse(UInt128(0x13579bdf02468ace8642fdb975310eca), bps)
    end
end
module FASTQ

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioCore: BioCore, isfilled
import BioSymbols
import BioSequences: BioSequences, ReaderIterTools
import TranscodingStreams: TranscodingStreams, TranscodingStream

include("quality.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")

end

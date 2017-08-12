module FASTQ

import Automa
import Automa.RegExp: @re_str
import BioCore.ReaderHelper: @pos, @mark, @unmark
import BioCore: BioCore, isfilled
import BioSequences
import BioSymbols
import CodecZlib
import TranscodingStreams: TranscodingStreams, TranscodingStream

include("quality.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")

end

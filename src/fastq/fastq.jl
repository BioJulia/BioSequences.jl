module FASTQ

import Automa
import Automa.RegExp: @re_str
import BioCore.ReaderHelper: @pos, @mark, @unmark, RecordIterator, RecordIteratorState, readrecord!
import BioCore: BioCore, isfilled
import BioSequences
import BioSymbols
import CodecZlib
import TranscodingStreams: TranscodingStream, NoopStream

include("quality.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")

end

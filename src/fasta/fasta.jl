# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import BioCore.Exceptions: missingerror
import BioCore.ReaderHelper: @pos, @mark, @unmark, RecordIterator, RecordIteratorState, readrecord!
import BioCore: BioCore, isfilled
import BioSequences
import CodecZlib
import TranscodingStreams: TranscodingStream, NoopStream

export description,
       identifier

include("record.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end

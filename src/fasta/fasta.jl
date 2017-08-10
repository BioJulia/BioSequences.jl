# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import BioCore: BioCore, isfilled
import BioCore.Exceptions: missingerror
import BioCore.ReaderHelper: @pos, @mark, @unmark
import BioSequences
import TranscodingStreams: TranscodingStreams, TranscodingStream
import CodecZlib

export description,
       identifier

include("record.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end

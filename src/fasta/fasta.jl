# FASTA File Format
# =================

module FASTA

import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos
import BioCore: BioCore, isfilled
import BioCore.Exceptions: missingerror
import BioSequences
import TranscodingStreams: TranscodingStreams, TranscodingStream
#import BufferedStreams

export description,
       identifier

include("record.jl")
include("index.jl")
include("reader.jl")
include("writer.jl")

end

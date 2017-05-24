module FASTQ

import Automa
import Automa.RegExp: @re_str
import BioCore: BioCore, isfilled
import BioSymbols
import BioSequences
import BufferedStreams
import BufferedStreams: BufferedInputStream

include("quality.jl")
include("record.jl")
include("reader.jl")
include("writer.jl")

end

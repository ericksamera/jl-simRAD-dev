using FASTX

struct Enzyme
    name::String
    pattern::String
    regex_pattern::Regex
    cut_site::Integer
    catalyze::Function
end

_DEGEN_NUC_TO_REGEX = Base.ImmutableDict(
    "W"=>"[A|T]",
    "S"=>"[C|G]",
    "M"=>"[A|C]",
    "K"=>"[G|T]",
    "R"=>"[A|G]",
    "Y"=>"[C|T]",
    "B"=>"[C|G|T]",
    "D"=>"[A|G|T]",
    "H"=>"[A|C|T]",
    "V"=>"[A|C|G]",
    "N"=>"[A|C|T|G]")

function new_Enzyme(name::String, pattern::String)
    regex_pattern = pattern
    for (key, value) in _DEGEN_NUC_TO_REGEX
        regex_pattern = replace(regex_pattern, key => value)
    end

    cut_site = findfirst("^", regex_pattern)[1]-2

    regex_pattern = replace(regex_pattern, "^"=>"")
    regex_pattern = replace(regex_pattern, "_"=>"")
    
    function catalyze(sequence::String)
        match_positions = vcat([0], [x.offset + cut_site for x in eachmatch(Regex(regex_pattern), sequence)], length(sequence))
        fragments = [(sequence[match_positions[i]+1:match_positions[i+1]]) for (i, value) in enumerate(1:length(match_positions)-1)]
        return fragments
    end

    return Enzyme(name, pattern, Regex(regex_pattern, "i"), cut_site, catalyze)
end

EcoRI = new_Enzyme("EcoRI", "G^AATT_C")

full_genome = "GCA_940337035.1_PGI_AGRIOTES_LIN_V1_genomic.fna"
test_seq = "16S.fasta"

open(FASTA.Reader, test_seq,) do reader
    for record in reader
        slices = EcoRI.catalyze(String(FASTX.sequence(record)))
        for slice in slices
            println(length(slice))
            println(slice)
        end
    end
end
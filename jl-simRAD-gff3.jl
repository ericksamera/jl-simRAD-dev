using FASTX
using ArgParse
using Printf

include("restriction_enzymes.jl")

struct FragmentInfo
    genome_size::Int
    bases_covered::Int
    num_fragments::Int
    total_fragments_before_filtering::Int
    fragments::Vector{RestrictionFragment}
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--min-size", "-m"
            help = "minimize size of fragments"
            arg_type = Int
            default = 300
        "--max-size", "-M"
            help = "minimize size of fragments"
            arg_type = Int
            default = 600

        "--pretty"
            help = "minimize size of fragments"
            arg_type = Bool
            default = true
        "--csv"
            help = "minimize size of fragments"
            arg_type = Bool
            default = true

        "--gff3"
            help = "minimize size of fragments"
            arg_type = Bool
            default = false

        "genome_path"
            help = "a positional argument"
            required = true
        "enzyme_1"
            help = "a positional argument"
            required = true
        "enzyme_2"
            help = "a positional argument"
            required = true
    end

    return parse_args(s)
end

function digest_sequence(sequence::String, enzyme1::Enzyme, enzyme2::Enzyme, min_size::Int, max_size::Int, chromosome_id::String)
    initial_fragments = [RestrictionFragment(sequence, "original", "original", 1, length(sequence), chromosome_id)]
    
    fragments_after_first_digest = catalyze(enzyme1, initial_fragments)

    bases_covered = 0
    total_fragments = 0
    num_fragments_after_filter = 0

    fragments_collected = []

    for fragment in fragments_after_first_digest
        fragments_after_second_digest = catalyze(enzyme2, [fragment])
        
        for fragment2 in fragments_after_second_digest
            if fragment2.left_end != fragment2.right_end
                len = length(fragment2.sequence)
                total_fragments += 1
                if min_size < len < max_size
                    bases_covered += len
                    num_fragments_after_filter += 1
                    push!(fragments_collected, fragment2)
                end
            end
        end
    end
    
    return bases_covered, num_fragments_after_filter, total_fragments, fragments_collected
end

function main()

    parsed_args = parse_commandline()

    if !haskey(enzyme_dict, parsed_args["enzyme_1"])
        println("Enzyme 1 is not recognized!")
    end

    if !haskey(enzyme_dict, parsed_args["enzyme_2"])
        println("Enzyme 2 is not recognized!")
    end

    println(stderr, "## ", parsed_args["enzyme_1"], "-",parsed_args["enzyme_2"])


    collected_fragments = []

    info = open(FASTA.Reader, parsed_args["genome_path"]) do reader
        total_genome_size = 0
        total_bases_covered = 0
        total_num_fragments = 0
        total_fragments_before_filtering = 0
        
        collected_fragments = []
        
        for record in reader
            sequence = String(FASTX.sequence(record))
            chromosome_id = FASTX.description(record)  # Or however you extract the chromosome ID
            total_genome_size += length(sequence)
    
            bases_covered, num_fragments, total_fragments, fragments = digest_sequence(
                sequence, 
                enzyme_dict[parsed_args["enzyme_1"]], 
                enzyme_dict[parsed_args["enzyme_2"]], 
                parsed_args["min-size"], 
                parsed_args["max-size"],
                String(FASTX.identifier(record))
            )

            total_bases_covered += bases_covered
            total_num_fragments += num_fragments
            total_fragments_before_filtering += total_fragments
            append!(collected_fragments, fragments)
        end
        
        FragmentInfo(total_genome_size, total_bases_covered, total_num_fragments, total_fragments_before_filtering, collected_fragments)
    end

    if parsed_args["csv"]
        @printf("%s-%s,%s,%s,%s,%s\n", parsed_args["enzyme_1"], parsed_args["enzyme_2"], info.bases_covered, info.num_fragments, info.total_fragments_before_filtering, info.bases_covered / info.genome_size * 100)
    end
    if parsed_args["pretty"]
        @printf("# %s-%s\n", parsed_args["enzyme_1"], parsed_args["enzyme_2"])
        println("# Genome size: ", info.genome_size)
        println("# Bases covered: ", info.bases_covered)
        println("# Number of fragments after filtering: ", info.num_fragments)
        println("# Number of fragments before filtering: ", info.total_fragments_before_filtering)
        println("# Percentage of bases covered: ", info.bases_covered / info.genome_size * 100, "%")
    end

    if parsed_args["gff3"]
        println("##gff-version 3")
        for fragment in info.fragments
            @printf("%s\t%s\t%s\t%d\t%d\t.\t.\t.\tID=fragment_%d;Name=fragment_%d\n", 
                    fragment.id, # Using the chromosome or contig ID from the fragment
                    "simRAD", 
                    "restriction_fragment",
                    fragment.start_pos,
                    fragment.end_pos,
                    fragment.start_pos, # Just using start_pos as a unique ID here
                    fragment.start_pos) 
        end
    end
end

main()
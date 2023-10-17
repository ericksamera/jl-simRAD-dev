using FASTX
using Combinatorics

include("restriction_enzymes.jl")

full_genome = "GCA_940337035.1_PGI_AGRIOTES_LIN_V1_genomic.fna"
test_seq = "16S.fasta"

size_filter = [300, 600]

for combo in combinations([Enzyme.name for Enzyme in restriction_enzymes_list], 2)
    println(join(combo, "-"))
end

# genome_size, bases_covered = open(FASTA.Reader, full_genome) do reader
#     local_genome_size = 0
#     local_bases_covered = 0

#     for record in reader
#         local_genome_size += length(FASTX.sequence(record))

#         initial_fragment = RestrictionFragment(String(FASTX.sequence(record)), "original", "original")
#         digest_1 = catalyze(EcoRI, [initial_fragment])
#         for digest_1_fragments in digest_1
#             digest_2 = catalyze(NcoI, [digest_1_fragments])
#             for digest_2_fragment in digest_2
#                 if digest_2_fragment.left_end == digest_2_fragment.right_end
#                     continue
#                 end

#                 if !(size_filter[1] < length(digest_2_fragment.sequence) < size_filter[2])
#                     continue
#                 end
                
#                 local_bases_covered += length(digest_2_fragment.sequence)
#             end
#         end
#     end

#     return local_genome_size, local_bases_covered
# end

# println(genome_size)
# println(bases_covered)
# println(bases_covered/genome_size * 100)
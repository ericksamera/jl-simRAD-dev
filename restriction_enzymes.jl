struct Enzyme
    name::String
    pattern::String
    regex_pattern::Regex
    cut_site::Int
end

struct RestrictionFragment
    sequence::String
    left_end::String
    right_end::String
end

const _DEGEN_NUC_TO_REGEX = Base.ImmutableDict(
    'A'=>"[A]",
    'T'=>"[T]",
    'C'=>"[C]",
    'G'=>"[G]",
    'W'=>"[AT]",
    'S'=>"[CG]",
    'M'=>"[AC]",
    'K'=>"[GT]",
    'R'=>"[AG]",
    'Y'=>"[CT]",
    'B'=>"[CGT]",
    'D'=>"[AGT]",
    'H'=>"[ACT]",
    'V'=>"[ACG]",
    'N'=>"[ACTG]")

function convert_pattern(pattern::String)
    regex_pattern = pattern
    for (key, value) in _DEGEN_NUC_TO_REGEX
        regex_pattern = replace(regex_pattern, key => value)
    end
    regex_pattern = replace(regex_pattern, "^"=>"", "_"=>"")
    return Regex(regex_pattern, "i")
end

function catalyze(enzyme::Enzyme, input_fragments::Vector{RestrictionFragment})
    output_fragments = RestrictionFragment[]

    for fragment in input_fragments
        match_positions = [0; [x.offset + enzyme.cut_site for x in eachmatch(enzyme.regex_pattern, fragment.sequence)]; length(fragment.sequence)]
        for i in 1:length(match_positions)-1
            sequence = fragment.sequence[match_positions[i]+1:match_positions[i+1]]
            left_end = i == 1 ? fragment.left_end : enzyme.name
            right_end = i == length(match_positions)-1 ? fragment.right_end : enzyme.name
            push!(output_fragments, RestrictionFragment(sequence, left_end, right_end))
        end
    end

    return output_fragments
end

function new_Enzyme(name::String, pattern::String)
    regex_pattern = convert_pattern(pattern)
    cut_site = findfirst("^", pattern).start - 2
    return Enzyme(name, pattern, regex_pattern, cut_site)
end


AatII=new_Enzyme("AatII", "G_ACGT^C")
#AbaSI=new_Enzyme("AbaSI", "hCNNNNNNNNNN_NN^")
Acc65I=new_Enzyme("Acc65I", "G^GTAC_C")
AccI=new_Enzyme("AccI", "GT^MK_AC")
AciI=new_Enzyme("AciI", "C^CG_C")
AclI=new_Enzyme("AclI", "AA^CG_TT")
AcuI=new_Enzyme("AcuI", "CTGAAGNNNNNNNNNNNNNNN_NN^")
AfeI=new_Enzyme("AfeI", "AGC_^GCT")
AflII=new_Enzyme("AflII", "C^TTAA_G")
AflIII=new_Enzyme("AflIII", "A^CRYG_T")
AgeI=new_Enzyme("AgeI", "A^CCGG_T")
AgeI_HF=new_Enzyme("AgeI-HF", "A^CCGG_T")
AhdI=new_Enzyme("AhdI", "GACNN_N^NNGTC")
AleI_v2=new_Enzyme("AleI-v2", "CACNN_^NNGTG")
AluI=new_Enzyme("AluI", "AG_^CT")
AlwI=new_Enzyme("AlwI", "GGATCNNNN^N_")
AlwNI=new_Enzyme("AlwNI", "CAG_NNN^CTG")
ApaI=new_Enzyme("ApaI", "G_GGCC^C")
ApaLI=new_Enzyme("ApaLI", "G^TGCA_C")
ApeKI=new_Enzyme("ApeKI", "G^CWG_C")
ApoI=new_Enzyme("ApoI", "R^AATT_Y")
ApoI_HF=new_Enzyme("ApoI-HF", "R^AATT_Y")
AscI=new_Enzyme("AscI", "GG^CGCG_CC")
AseI=new_Enzyme("AseI", "AT^TA_AT")
AsiSI=new_Enzyme("AsiSI", "GCG_AT^CGC")
AvaI=new_Enzyme("AvaI", "C^YCGR_G")
AvaII=new_Enzyme("AvaII", "G^GWC_C")
AvrII=new_Enzyme("AvrII", "C^CTAG_G")
BaeGI=new_Enzyme("BaeGI", "G_KGCM^C")
BaeI=new_Enzyme("BaeI", "_NNNNN^NNNNNNNNNNNACNNNNGTAYCNNNNNNNN_NNNNN^")
BamHI=new_Enzyme("BamHI", "G^GATC_C")
BamHI_HF=new_Enzyme("BamHI-HF", "G^GATC_C")
BanI=new_Enzyme("BanI", "G^GYRC_C")
BanII=new_Enzyme("BanII", "G_RGCY^C")
BbsI=new_Enzyme("BbsI", "GAAGACNN^NNNN_")
BbsI_HF=new_Enzyme("BbsI-HF", "GAAGACNN^NNNN_")
BbvCI=new_Enzyme("BbvCI", "CC^TCA_GC")
BbvI=new_Enzyme("BbvI", "GCAGCNNNNNNNNN^NNNN_")
BccI=new_Enzyme("BccI", "CCATCNNNN^N_")
BceAI=new_Enzyme("BceAI", "ACGGCNNNNNNNNNNNNN^NN_")
BcgI=new_Enzyme("BcgI", "_NN^NNNNNNNNNNNCGANNNNNNTGCN_NN^")
BciVI=new_Enzyme("BciVI", "GTATCCNNNNNN_N^")
BclI=new_Enzyme("BclI", "T^GATC_A")
BclI_HF=new_Enzyme("BclI-HF", "T^GATC_A")
BcoDI=new_Enzyme("BcoDI", "GTCTCN^NNNN_")
BfaI=new_Enzyme("BfaI", "C^TA_G")
BfuAI=new_Enzyme("BfuAI", "ACCTGCNNNN^NNNN_")
BglI=new_Enzyme("BglI", "GCCN_NNN^NGGC")
BglII=new_Enzyme("BglII", "A^GATC_T")
BlpI=new_Enzyme("BlpI", "GC^TNA_GC")
BmgBI=new_Enzyme("BmgBI", "CAC_^GTC")
BmrI=new_Enzyme("BmrI", "ACTGGGNNNN_N^")
BmtI=new_Enzyme("BmtI", "G_CTAG^C")
BmtI_HF=new_Enzyme("BmtI-HF", "G_CTAG^C")
BpmI=new_Enzyme("BpmI", "CTGGAGNNNNNNNNNNNNNNN_NN^")
Bpu10I=new_Enzyme("Bpu10I", "CC^TNA_GC")
BpuEI=new_Enzyme("BpuEI", "CTTGAGNNNNNNNNNNNNNNN_NN^")
BsaAI=new_Enzyme("BsaAI", "YAC_^GTR")
BsaBI=new_Enzyme("BsaBI", "GATNN_^NNATC")
BsaHI=new_Enzyme("BsaHI", "GR^CG_YC")
BsaI_HFv2=new_Enzyme("BsaI-HFv2", "GGTCTCN^NNNN_")
BsaJI=new_Enzyme("BsaJI", "C^CNNG_G")
BsaWI=new_Enzyme("BsaWI", "W^CCGG_W")
BsaXI=new_Enzyme("BsaXI", "_NNN^NNNNNNNNNNACNNNNNCTCCNNNNNNNN_NNN^")
BseRI=new_Enzyme("BseRI", "GAGGAGNNNNNNNNN_NN^")
BseYI=new_Enzyme("BseYI", "C^CCAG_C")
BsgI=new_Enzyme("BsgI", "GTGCAGNNNNNNNNNNNNNNN_NN^")
BsiEI=new_Enzyme("BsiEI", "CG_RY^CG")
BsiHKAI=new_Enzyme("BsiHKAI", "G_WGCW^C")
BsiWI=new_Enzyme("BsiWI", "C^GTAC_G")
BsiWI_HF=new_Enzyme("BsiWI-HF", "C^GTAC_G")
BslI=new_Enzyme("BslI", "CCNN_NNN^NNGG")
BsmAI=new_Enzyme("BsmAI", "GTCTCN^NNNN_")
BsmBI_v2=new_Enzyme("BsmBI-v2", "CGTCTCN^NNNN_")
BsmFI=new_Enzyme("BsmFI", "GGGACNNNNNNNNNNN^NNNN_")
BsmI=new_Enzyme("BsmI", "GAATG_CN^")
BsoBI=new_Enzyme("BsoBI", "C^YCGR_G")
Bsp1286I=new_Enzyme("Bsp1286I", "G_DGCH^C")
BspCNI=new_Enzyme("BspCNI", "CTCAGNNNNNNNN_NN^")
BspDI=new_Enzyme("BspDI", "AT^CG_AT")
BspEI=new_Enzyme("BspEI", "T^CCGG_A")
BspHI=new_Enzyme("BspHI", "T^CATG_A")
BspMI=new_Enzyme("BspMI", "ACCTGCNNNN^NNNN_")
BspQI=new_Enzyme("BspQI", "GCTCTTCN^NNN_")
BsrBI=new_Enzyme("BsrBI", "CCG_^CTC")
BsrDI=new_Enzyme("BsrDI", "GCAATG_NN^")
BsrFI_v2=new_Enzyme("BsrFI-v2", "R^CCGG_Y")
BsrGI=new_Enzyme("BsrGI", "T^GTAC_A")
BsrGI_HF=new_Enzyme("BsrGI-HF", "T^GTAC_A")
BsrI=new_Enzyme("BsrI", "ACTG_GN^")
BssHII=new_Enzyme("BssHII", "G^CGCG_C")
BssSI_v2=new_Enzyme("BssSI-v2", "C^ACGA_G")
BstAPI=new_Enzyme("BstAPI", "GCAN_NNN^NTGC")
BstBI=new_Enzyme("BstBI", "TT^CG_AA")
BstEII=new_Enzyme("BstEII", "G^GTNAC_C")
BstEII_HF=new_Enzyme("BstEII-HF", "G^GTNAC_C")
BstNI=new_Enzyme("BstNI", "CC^W_GG")
BstUI=new_Enzyme("BstUI", "CG_^CG")
BstXI=new_Enzyme("BstXI", "CCAN_NNNN^NTGG")
BstYI=new_Enzyme("BstYI", "R^GATC_Y")
BstZ17I_HF=new_Enzyme("BstZ17I-HF", "GTA_^TAC")
Bsu36I=new_Enzyme("Bsu36I", "CC^TNA_GG")
BtgI=new_Enzyme("BtgI", "C^CRYG_G")
BtgZI=new_Enzyme("BtgZI", "GCGATGNNNNNNNNNNN^NNNN_")
BtsCI=new_Enzyme("BtsCI", "GGATG_NN^")
BtsI_v2=new_Enzyme("BtsI-v2", "GCAGTG_NN^")
BtsIMutI=new_Enzyme("BtsIMutI", "CAGTG_NN^")
Cac8I=new_Enzyme("Cac8I", "GCN_^NGC")
ClaI=new_Enzyme("ClaI", "AT^CG_AT")
CspCI=new_Enzyme("CspCI", "_NN^NNNNNNNNNNNNCAANNNNNGTGGNNNNNNNNNNN_NN^")
CviAII=new_Enzyme("CviAII", "C^AT_G")
CviKI_1=new_Enzyme("CviKI-1", "RG_^CY")
CviQI=new_Enzyme("CviQI", "G^TA_C")
DdeI=new_Enzyme("DdeI", "C^TNA_G")
DpnI=new_Enzyme("DpnI", "GA_^TC")
DpnII=new_Enzyme("DpnII", "^GATC_")
DraI=new_Enzyme("DraI", "TTT_^AAA")
DraIII_HF=new_Enzyme("DraIII-HF", "CAC_NNN^GTG")
DrdI=new_Enzyme("DrdI", "GACNN_NN^NNGTC")
EaeI=new_Enzyme("EaeI", "Y^GGCC_R")
EagI_HF=new_Enzyme("EagI-HF", "C^GGCC_G")
EarI=new_Enzyme("EarI", "CTCTTCN^NNN_")
EciI=new_Enzyme("EciI", "GGCGGANNNNNNNNNN_NN^")
Eco53kI=new_Enzyme("Eco53kI", "GAG_^CTC")
EcoNI=new_Enzyme("EcoNI", "CCTNN^N_NNAGG")
EcoO109I=new_Enzyme("EcoO109I", "RG^GNC_CY")
EcoRI=new_Enzyme("EcoRI", "G^AATT_C")
EcoRI_HF=new_Enzyme("EcoRI-HF", "G^AATT_C")
EcoRV=new_Enzyme("EcoRV", "GAT_^ATC")
EcoRV_HF=new_Enzyme("EcoRV-HF", "GAT_^ATC")
Esp3I=new_Enzyme("Esp3I", "CGTCTCN^NNNN_")
FatI=new_Enzyme("FatI", "^CATG_")
FauI=new_Enzyme("FauI", "CCCGCNNNN^NN_")
Fnu4HI=new_Enzyme("Fnu4HI", "GC^N_GC")
FokI=new_Enzyme("FokI", "GGATGNNNNNNNNNN^NNNN_")
FseI=new_Enzyme("FseI", "GG_CCGG^CC")
FspEI=new_Enzyme("FspEI", "CCNNNNNNNNNNNNN^NNNN_")
FspI=new_Enzyme("FspI", "TGC_^GCA")
HaeII=new_Enzyme("HaeII", "R_GCGC^Y")
HaeIII=new_Enzyme("HaeIII", "GG_^CC")
HgaI=new_Enzyme("HgaI", "GACGCNNNNNN^NNNNN_")
HhaI=new_Enzyme("HhaI", "G_CG^C")
HinP1I=new_Enzyme("HinP1I", "G^CG_C")
HincII=new_Enzyme("HincII", "GTY_^RAC")
HindIII=new_Enzyme("HindIII", "A^AGCT_T")
HindIII_HF=new_Enzyme("HindIII-HF", "A^AGCT_T")
HinfI=new_Enzyme("HinfI", "G^ANT_C")
HpaI=new_Enzyme("HpaI", "GTT_^AAC")
HpaII=new_Enzyme("HpaII", "C^CG_G")
HphI=new_Enzyme("HphI", "GGTGANNNNNNNN_N^")
Hpy166II=new_Enzyme("Hpy166II", "GTN_^NAC")
Hpy188I=new_Enzyme("Hpy188I", "TC_N^GA")
Hpy188III=new_Enzyme("Hpy188III", "TC^NN_GA")
Hpy99I=new_Enzyme("Hpy99I", "_CGWCG^")
HpyAV=new_Enzyme("HpyAV", "CCTTCNNNNNN_N^")
HpyCH4III=new_Enzyme("HpyCH4III", "AC_N^GT")
HpyCH4IV=new_Enzyme("HpyCH4IV", "A^CG_T")
HpyCH4V=new_Enzyme("HpyCH4V", "TG_^CA")
I_CeuI=new_Enzyme("I-CeuI", "TAACTATAACGGTC_CTAA^GGTAGCGAA")
I_SceI=new_Enzyme("I-SceI", "TAGGG_ATAA^CAGGGTAAT")
KasI=new_Enzyme("KasI", "G^GCGC_C")
KpnI=new_Enzyme("KpnI", "G_GTAC^C")
KpnI_HF=new_Enzyme("KpnI-HF", "G_GTAC^C")
LpnPI=new_Enzyme("LpnPI", "CCDGNNNNNNNNNNN^NNNN_")
MboI=new_Enzyme("MboI", "^GATC_")
MboII=new_Enzyme("MboII", "GAAGANNNNNNNN_N^")
MfeI=new_Enzyme("MfeI", "C^AATT_G")
MfeI_HF=new_Enzyme("MfeI-HF", "C^AATT_G")
MluCI=new_Enzyme("MluCI", "^AATT_")
MluI=new_Enzyme("MluI", "A^CGCG_T")
MluI_HF=new_Enzyme("MluI-HF", "A^CGCG_T")
MlyI=new_Enzyme("MlyI", "GAGTCNNNNNN_^")
MmeI=new_Enzyme("MmeI", "TCCRACNNNNNNNNNNNNNNNNNNN_NN^")
MnlI=new_Enzyme("MnlI", "CCTCNNNNNNN_N^")
MscI=new_Enzyme("MscI", "TGG_^CCA")
MseI=new_Enzyme("MseI", "T^TA_A")
MslI=new_Enzyme("MslI", "CAYNN_^NNRTG")
MspA1I=new_Enzyme("MspA1I", "CMG_^CKG")
MspI=new_Enzyme("MspI", "C^CG_G")
MspJI=new_Enzyme("MspJI", "CNNRNNNNNNNNNN^NNNN_")
MwoI=new_Enzyme("MwoI", "GCNN_NNN^NNGC")
NaeI=new_Enzyme("NaeI", "GCC_^GGC")
NarI=new_Enzyme("NarI", "GG^CG_CC")
NciI=new_Enzyme("NciI", "CC^S_GG")
NcoI=new_Enzyme("NcoI", "C^CATG_G")
NcoI_HF=new_Enzyme("NcoI-HF", "C^CATG_G")
NdeI=new_Enzyme("NdeI", "CA^TA_TG")
NgoMIV=new_Enzyme("NgoMIV", "G^CCGG_C")
NheI_HF=new_Enzyme("NheI-HF", "G^CTAG_C")
NlaIII=new_Enzyme("NlaIII", "_CATG^")
NlaIV=new_Enzyme("NlaIV", "GGN_^NCC")
NmeAIII=new_Enzyme("NmeAIII", "GCCGAGNNNNNNNNNNNNNNNNNNNN_N^")
NotI=new_Enzyme("NotI", "GC^GGCC_GC")
NotI_HF=new_Enzyme("NotI-HF", "GC^GGCC_GC")
NruI=new_Enzyme("NruI", "TCG_^CGA")
NruI_HF=new_Enzyme("NruI-HF", "TCG_^CGA")
NsiI=new_Enzyme("NsiI", "A_TGCA^T")
NsiI_HF=new_Enzyme("NsiI-HF", "A_TGCA^T")
NspI=new_Enzyme("NspI", "R_CATG^Y")
Nt_AlwI=new_Enzyme("Nt.AlwI", "GGATCNNNN^",)
Nt_BbvCI=new_Enzyme("Nt.BbvCI", "CC^TCAGC",)
Nt_BsmAI=new_Enzyme("Nt.BsmAI", "GTCTCN^",)
Nt_BspQI=new_Enzyme("Nt.BspQI", "GCTCTTCN^",)
Nt_BstNBI=new_Enzyme("Nt.BstNBI", "GAGTCNNNN^",)
Nt_CviPII=new_Enzyme("Nt.CviPII", "^CCD",)
PI_PspI=new_Enzyme("PI-PspI", "TGGCAAACAGCTA_TTAT^GGGTATTATGGGT")
PI_SceI=new_Enzyme("PI-SceI", "ATCTATGTCGG_GTGC^GGAGAAAGAGGTAATGAAATGG")
PacI=new_Enzyme("PacI", "TTA_AT^TAA")
PaeR7I=new_Enzyme("PaeR7I", "C^TCGA_G")
PaqCI=new_Enzyme("PaqCI", "CACCTGCNNNN^NNNN_")
PciI=new_Enzyme("PciI", "A^CATG_T")
PflFI=new_Enzyme("PflFI", "GACN^N_NGTC")
PflMI=new_Enzyme("PflMI", "CCAN_NNN^NTGG")
PleI=new_Enzyme("PleI", "GAGTCNNNN^N_")
PluTI=new_Enzyme("PluTI", "G_GCGC^C")
PmeI=new_Enzyme("PmeI", "GTTT_^AAAC")
PmlI=new_Enzyme("PmlI", "CAC_^GTG")
PpuMI=new_Enzyme("PpuMI", "RG^GWC_CY")
PshAI=new_Enzyme("PshAI", "GACNN_^NNGTC")
PsiI_v2=new_Enzyme("PsiI-v2", "TTA_^TAA")
PspGI=new_Enzyme("PspGI", "^CCWGG_")
PspOMI=new_Enzyme("PspOMI", "G^GGCC_C")
PspXI=new_Enzyme("PspXI", "VC^TCGA_GB")
PstI=new_Enzyme("PstI", "C_TGCA^G")
PstI_HF=new_Enzyme("PstI-HF", "C_TGCA^G")
PvuI=new_Enzyme("PvuI", "CG_AT^CG")
PvuI_HF=new_Enzyme("PvuI-HF", "CG_AT^CG")
PvuII=new_Enzyme("PvuII", "CAG_^CTG")
PvuII_HF=new_Enzyme("PvuII-HF", "CAG_^CTG")
RsaI=new_Enzyme("RsaI", "GT_^AC")
RsrII=new_Enzyme("RsrII", "CG^GWC_CG")
SacI=new_Enzyme("SacI", "G_AGCT^C")
SacI_HF=new_Enzyme("SacI-HF", "G_AGCT^C")
SacII=new_Enzyme("SacII", "CC_GC^GG")
SalI=new_Enzyme("SalI", "G^TCGA_C")
SalI_HF=new_Enzyme("SalI-HF", "G^TCGA_C")
SapI=new_Enzyme("SapI", "GCTCTTCN^NNN_")
Sau3AI=new_Enzyme("Sau3AI", "^GATC_")
Sau96I=new_Enzyme("Sau96I", "G^GNC_C")
SbfI=new_Enzyme("SbfI", "CC_TGCA^GG")
SbfI_HF=new_Enzyme("SbfI-HF", "CC_TGCA^GG")
ScaI_HF=new_Enzyme("ScaI-HF", "AGT_^ACT")
ScrFI=new_Enzyme("ScrFI", "CC^N_GG")
SexAI=new_Enzyme("SexAI", "A^CCWGG_T")
SfaNI=new_Enzyme("SfaNI", "GCATCNNNNNN^NNNN_")
SfcI=new_Enzyme("SfcI", "C^TRYA_G")
SfiI=new_Enzyme("SfiI", "GGCCN_NNN^NGGCC")
SfoI=new_Enzyme("SfoI", "GGC_^GCC")
SgrAI=new_Enzyme("SgrAI", "CR^CCGG_YG")
SmaI=new_Enzyme("SmaI", "CCC_^GGG")
SmlI=new_Enzyme("SmlI", "C^TYRA_G")
SnaBI=new_Enzyme("SnaBI", "TAC_^GTA")
SpeI=new_Enzyme("SpeI", "A^CTAG_T")
SpeI_HF=new_Enzyme("SpeI-HF", "A^CTAG_T")
SphI=new_Enzyme("SphI", "G_CATG^C")
SphI_HF=new_Enzyme("SphI-HF", "G_CATG^C")
SrfI=new_Enzyme("SrfI", "GCCC_^GGGC")
SspI=new_Enzyme("SspI", "AAT_^ATT")
SspI_HF=new_Enzyme("SspI-HF", "AAT_^ATT")
StuI=new_Enzyme("StuI", "AGG_^CCT")
StyD4I=new_Enzyme("StyD4I", "^CCNGG_")
StyI_HF=new_Enzyme("StyI-HF", "C^CWWG_G")
SwaI=new_Enzyme("SwaI", "ATTT_^AAAT")
TaqI_v2=new_Enzyme("TaqI-v2", "T^CG_A")
TfiI=new_Enzyme("TfiI", "G^AWT_C")
TseI=new_Enzyme("TseI", "G^CWG_C")
Tsp45I=new_Enzyme("Tsp45I", "^GTSAC_")
TspMI=new_Enzyme("TspMI", "C^CCGG_G")
TspRI=new_Enzyme("TspRI", "_NNCASTGNN^")
Tth111I=new_Enzyme("Tth111I", "GACN^N_NGTC")
XbaI=new_Enzyme("XbaI", "T^CTAG_A")
XcmI=new_Enzyme("XcmI", "CCANNNN_N^NNNNTGG")
XhoI=new_Enzyme("XhoI", "C^TCGA_G")
XmaI=new_Enzyme("XmaI", "C^CCGG_G")
XmnI=new_Enzyme("XmnI", "GAANN_^NNTTC")
ZraI=new_Enzyme("ZraI", "GAC_^GTC")

restriction_enzymes_list = [AatII, Acc65I, AccI, AciI, AclI, AcuI, AfeI, AflII, AflIII, AgeI, AgeI_HF, AhdI, AleI_v2, AluI, AlwI, AlwNI, ApaI, ApaLI, ApeKI, ApoI, ApoI_HF, AscI, AseI, AsiSI, AvaI, AvaII, AvrII, BaeGI, BaeI, BamHI, BamHI_HF, BanI, BanII, BbsI, BbsI_HF, BbvCI, BbvI, BccI, BceAI, BcgI, BciVI, BclI, BclI_HF, BcoDI, BfaI, BfuAI, BglI, BglII, BlpI, BmgBI, BmrI, BmtI, BmtI_HF, BpmI, Bpu10I, BpuEI, BsaAI, BsaBI, BsaHI, BsaI_HFv2, BsaJI, BsaWI, BsaXI, BseRI, BseYI, BsgI, BsiEI, BsiHKAI, BsiWI, BsiWI_HF, BslI, BsmAI, BsmBI_v2, BsmFI, BsmI, BsoBI, Bsp1286I, BspCNI, BspDI, BspEI, BspHI, BspMI, BspQI, BsrBI, BsrDI, BsrFI_v2, BsrGI, BsrGI_HF, BsrI, BssHII, BssSI_v2, BstAPI, BstBI, BstEII, BstEII_HF, BstNI, BstUI, BstXI, BstYI, BstZ17I_HF, Bsu36I, BtgI, BtgZI, BtsCI, BtsI_v2, BtsIMutI, Cac8I, ClaI, CspCI, CviAII, CviKI_1, CviQI, DdeI, DpnI, DpnII, DraI, DraIII_HF, DrdI, EaeI, EagI_HF, EarI, EciI, Eco53kI, EcoNI, EcoO109I, EcoRI, EcoRI_HF, EcoRV, EcoRV_HF, Esp3I, FatI, FauI, Fnu4HI, FokI, FseI, FspEI, FspI, HaeII, HaeIII, HgaI, HhaI, HinP1I, HincII, HindIII, HindIII_HF, HinfI, HpaI, HpaII, HphI, Hpy166II, Hpy188I, Hpy188III, Hpy99I, HpyAV, HpyCH4III, HpyCH4IV, HpyCH4V, I_CeuI, I_SceI, KasI, KpnI, KpnI_HF, LpnPI, MboI, MboII, MfeI, MfeI_HF, MluCI, MluI, MluI_HF, MlyI, MmeI, MnlI, MscI, MseI, MslI, MspA1I, MspI, MspJI, MwoI, NaeI, NarI, NciI, NcoI, NcoI_HF, NdeI, NgoMIV, NheI_HF, NlaIII, NlaIV, NmeAIII, NotI, NotI_HF, NruI, NruI_HF, NsiI, NsiI_HF, NspI, Nt_AlwI, Nt_BbvCI, Nt_BsmAI, Nt_BspQI, Nt_BstNBI, Nt_CviPII, PI_PspI, PI_SceI, PacI, PaeR7I, PaqCI, PciI, PflFI, PflMI, PleI, PluTI, PmeI, PmlI, PpuMI, PshAI, PsiI_v2, PspGI, PspOMI, PspXI, PstI, PstI_HF, PvuI, PvuI_HF, PvuII, PvuII_HF, RsaI, RsrII, SacI, SacI_HF, SacII, SalI, SalI_HF, SapI, Sau3AI, Sau96I, SbfI, SbfI_HF, ScaI_HF, ScrFI, SexAI, SfaNI, SfcI, SfiI, SfoI, SgrAI, SmaI, SmlI, SnaBI, SpeI, SpeI_HF, SphI, SphI_HF, SrfI, SspI, SspI_HF, StuI, StyD4I, StyI_HF, SwaI, TaqI_v2, TfiI, TseI, Tsp45I, TspMI, TspRI, Tth111I, XbaI, XcmI, XhoI, XmaI, XmnI, ZraI]
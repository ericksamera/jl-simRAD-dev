from Bio import SeqIO

__description__ =\
"""
Purpose: Module for creating a streamlined restriction enzyme.
"""
__author__ = "Erick Samera"
__version__ = "1.0.0"
__comments__ = "fast enough"
# --------------------------------------------------
import re
# --------------------------------------------------
class RestrictionEnzyme:
    def __init__(self, name: str, cut_sequence: str, overhang: str, incubation_temp: int, buffer: str, time_saver: bool, heat_inactivated: bool) -> None:
        """Create a restriction enzyme"""
        _degenerate_nucleotides = {
            'W': '[A|T]',
            'S': '[C|G]',
            'M': '[A|C]',
            'K': '[G|T]',
            'R': '[A|G]',
            'Y': '[C|T]',
            'B': '[C|G|T]',
            'D': '[A|G|T]',
            'H': '[A|C|T]',
            'V': '[A|C|G]',
            'N': '[A|C|T|G]'
        }

        self.name: str = str(name)
        self.overhang: str = str(overhang)
        self.incubation_temp: int = int(incubation_temp)
        self.buffer: str = str(buffer)
        self.time_saver: bool = True if time_saver == 'True' else False
        self.heat_inactivated: bool = True if heat_inactivated == 'True' else False

        self.regex_restriction_site = cut_sequence
        for degenenerate_nucleotide, regex_equiv in _degenerate_nucleotides.items():
            self.regex_restriction_site = self.regex_restriction_site.replace(degenenerate_nucleotide, regex_equiv)
        self.regex_restriction_site = self.regex_restriction_site.replace('^', '').replace('_', '')

        self.cut_sequence_str: str = cut_sequence
        self.recognition_site: str = cut_sequence.replace('^', '').replace('_', '')
        self.cut_site: int = cut_sequence.index('^')
        return None
    def __repr__(self):
        return f"<{self.name}: {self.cut_sequence_str}, {self.regex_restriction_site}>"
    def catalyze(self, target_sequence: str):
        """Perform catalysis on a sequence and return the fragments"""
        re_pattern = re.compile(self.regex_restriction_site, flags=re.IGNORECASE)
        re_matches = re.finditer(re_pattern, str(target_sequence))
        match_positions: list = [0] + [m.start()+self.cut_site for m in re_matches] + [len(str(target_sequence))]
        fragments: list = [(target_sequence[match_positions[i]:match_positions[i+1]]) for i, _ in enumerate(match_positions) if i < len(match_positions)-1]
        return fragments
    def fast_catalyze(self, target_sequence: str):
        """Perform catalysis on a sequence and return the fragments"""
        re_pattern = re.compile(self.regex_restriction_site, flags=re.IGNORECASE)
        re_matches = re.finditer(re_pattern, str(target_sequence))
        match_positions: list = [0] + [m.start()+self.cut_site for m in re_matches] + [len(str(target_sequence))]
        fragments: list = [(match_positions[i], match_positions[i+1]) for i, _ in enumerate(match_positions) if i < len(match_positions)-1]
        return fragments

EcoRI=RestrictionEnzyme('EcoRI', 'G^AATT_C', '5′', 37, 'NEB U', True, True)
reader = SeqIO.parse('GCA_940337035.1_PGI_AGRIOTES_LIN_V1_genomic.fna', 'fasta')

for record in reader:
    EcoRI.catalyze(record.seq)
    pass
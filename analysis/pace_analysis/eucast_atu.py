"""EUCAST areas of technical uncertainty (ATU) for A. fumigatus. The rule is defined only here.

Source: EUCAST antifungal breakpoint table v12.1 (2026-04-10); its Aspergillus table is identical to v12.0
(2025-06-26). For A. fumigatus, two MICs are areas of technical uncertainty:
  - isavuconazole 2 mg/L (note 1): report S if voriconazole is wild type (<= 1 mg/L), otherwise R and refer for CYP51A
    sequencing;
  - posaconazole 0.25 mg/L (note 5): report S if itraconazole is S (<= 1 mg/L), otherwise R and refer.
The ATU is not a category. The yeast table has no ATU.
pace_tools/dictionaries.py has no ATU (S and R wells are adjacent), so these MICs count as R during threshold
optimisation; the paper's categories, error counts and figures apply this rule afterwards.
"""

ATU_MARK = '‡'

# (species code, plate row) -> rule. 'decided_by' is the plate row of the drug whose MIC decides the category.
EUCAST_ATU = {
    ('ASFU', 'E'): {'name': 'isavuconazole', 'mic': 2,
                    'decided_by': 'B', 'decider_name': 'voriconazole', 's_if_le': 1},
    ('ASFU', 'C'): {'name': 'posaconazole', 'mic': 0.25,
                    'decided_by': 'D', 'decider_name': 'itraconazole', 's_if_le': 1},
}


def resolve(rule, mic, decider_mic):
    """Return 'S' or 'R' for an MIC in the ATU of `rule`, or None if `mic` is not the ATU value.

    mic, decider_mic: numeric MICs in mg/L. An off-scale MIC is passed as one dilution beyond the tested range.
    decider_mic is None when the deciding drug has no MIC; that raises, because the ATU cannot be resolved.
    """
    if mic is None or mic != rule['mic']:
        return None
    if decider_mic is None:
        raise RuntimeError(f"{rule['name']} {rule['mic']:g} mg/L is in the EUCAST ATU, but there is no "
                           f"{rule['decider_name']} MIC to resolve it")
    return 'S' if decider_mic <= rule['s_if_le'] else 'R'

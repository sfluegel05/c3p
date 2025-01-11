"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or structurally related compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sufficient number of chlorine or bromine atoms
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 3:  # Lowered to capture more true positives
        return False, "Insufficient halogen atoms for a likely related compound."

    # SMARTS pattern for polychlorinated biphenyls (PCBs)
    biphenyl_pattern = Chem.MolFromSmarts('c1cc(c(cc1)-c2ccccc2)[Cl,Br]')
    
    # SMARTS patterns for polychlorinated dibenzodioxins (PCDDs)
    dioxin_patterns = [
        Chem.MolFromSmarts('Oc1cc2Oc3ccccc3Oc2cc1'),  # Allow variability in Cl positions
        Chem.MolFromSmarts('c1cc2Oc3ccc(Cl)cc3cc2cc1')
    ]
    
    # SMARTS patterns for polychlorinated dibenzofurans (PCDFs)
    furan_patterns = [
        Chem.MolFromSmarts('c1cc2oc3cc(Cl)ccc3c2cc1'),
        Chem.MolFromSmarts('c1cc2oc3ccc(Cl)c2cc3cc1')
    ]
    
    # Check for presence of these patterns
    matches_dioxin_furan = any(mol.HasSubstructMatch(pattern) for pattern in dioxin_patterns + furan_patterns)
    matches_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    # Exclude common false patterns
    common_false_patterns = [
        Chem.MolFromSmarts('[OH]'),  # Hydroxy groups are not typical in PCBs, PCDD/Fs
    ]
    false_positive_match = any(mol.HasSubstructMatch(pattern) for pattern in common_false_patterns)

    if (matches_dioxin_furan or matches_biphenyl) and not false_positive_match:
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern."

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern."
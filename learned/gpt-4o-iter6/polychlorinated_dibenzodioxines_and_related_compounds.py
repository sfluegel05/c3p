"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or structurally related compound such as PCBs, PCDDs, or PCDFs based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Check for sufficient number of chlorine or bromine atoms
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 3:
        return (False, "Insufficient halogen atoms for a likely related compound.")

    # Common substructure patterns
    biphenyl_pattern = Chem.MolFromSmarts('C1=CC=C(C=C1)-C1=CC=CC=C1')
    dioxin_pattern = Chem.MolFromSmarts('Oc1c2ccccc2Oc2ccccc12')
    furan_pattern = Chem.MolFromSmarts('c1cc2oc3ccccc3c2cc1')
    
    # Check against pattern matches with oxygen bridging as optional
    matches_dioxin_furan = any(mol.HasSubstructMatch(pattern) for pattern in [dioxin_pattern, furan_pattern])
    matches_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    # If matches to any common false-pattern, exclude these
    false_patterns = [Chem.MolFromSmarts('[OH]'), Chem.MolFromSmarts('C1=CC=C(C=C1)N')]
    false_positive_match = any(mol.HasSubstructMatch(pattern) for pattern in false_patterns)

    if (matches_dioxin_furan or matches_biphenyl) and not false_positive_match:
        return (True, "Matches polychlorinated dibenzodioxin or related compound pattern.")

    return (False, "Does not match polychlorinated dibenzodioxin or related compound pattern.")
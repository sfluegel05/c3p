"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or a related compound based on its SMILES string.
    
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

    # Check for sufficient halogen atoms (Cl or Br)
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 3:
        return False, "Not enough halogen atoms for a polychlorinated or polybrominated compound"

    # SMARTS patterns for checking the relevant structures
    pcb_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')  # Biphenyl structure
    dioxin_pattern = Chem.MolFromSmarts('c1cc2Oc3cc(Cl)ccc3Oc2cc1')  # Chlorinated dibenzodioxin
    furan_pattern = Chem.MolFromSmarts('c1cc2oc3cc(Cl)ccc3c2cc1')   # Chlorinated dibenzofuran

    # SMARTS pattern for a more general polychlorinated aromatic compound
    polychlorinated_aromatics_pattern = Chem.MolFromSmarts('c1cc(-c2ccccc2)cc1')

    # Check for presence of these patterns in the molecule
    if (mol.HasSubstructMatch(pcb_pattern) or
        mol.HasSubstructMatch(dioxin_pattern) or
        mol.HasSubstructMatch(furan_pattern) or
        mol.HasSubstructMatch(polychlorinated_aromatics_pattern)):
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern"

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern"

# Example to test the function
smiles_example = "Clc1cc2oc3ccccc3c2cc1"
result, reason = is_polychlorinated_dibenzodioxines_and_related_compounds(smiles_example)
print(f"Result: {result}, Reason: {reason}")
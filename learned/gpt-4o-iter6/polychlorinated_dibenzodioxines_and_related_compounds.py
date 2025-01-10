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
    if halogen_count < 2:
        return False, "Not enough halogen atoms for a relevant compound"

    # Define SMARTS pattern for biphenyl structure potentially with various chlorination
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')
    
    # Define SMARTS pattern for common polychlorinated dibenzodioxins
    dioxin_patterns = [
        Chem.MolFromSmarts('c1cc2Oc3ccc(cc3Oc2cc1)Cl'),  # Generic dioxin pattern
        Chem.MolFromSmarts('c1cc2Oc3cc(Cl)c(Cl)cc3Oc2cc1')  # Specific chlorination
    ]
    
    # Define SMARTS pattern for common polychlorinated dibenzofurans
    furan_patterns = [
        Chem.MolFromSmarts('c1cc2oc3ccc(cc3)c2cc1'),
        Chem.MolFromSmarts('c1cc2oc3cc(Cl)c(Cl)cc3c2cc1')
    ]

    # Check for presence of these patterns in the molecule
    if (mol.HasSubstructMatch(biphenyl_pattern) or
        any(mol.HasSubstructMatch(pattern) for pattern in dioxin_patterns) or
        any(mol.HasSubstructMatch(pattern) for pattern in furan_patterns)):
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern"

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern"
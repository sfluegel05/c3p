"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue
"""

from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    They typically have a core structure similar to canonical nucleobases (adenine, guanine, cytosine, thymine, uracil)
    with possible modifications, but without attached sugar or phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import AllChem

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for canonical nucleobases
    nucleobase_patterns = {
        'adenine': Chem.MolFromSmarts('c1nc2[nH]cnc-2nc1'),
        'guanine': Chem.MolFromSmarts('c1[nH]c2c(n1)nc(=O)[nH]c2'),
        'cytosine': Chem.MolFromSmarts('Nc1nc[nH]cc1=O'),
        'thymine': Chem.MolFromSmarts('O=C1NC(=O)C=C[NH]1'),
        'uracil': Chem.MolFromSmarts('O=C1NC(=O)C=CC1=O')
    }

    # Define SMARTS patterns for sugar (ribose/deoxyribose)
    sugar_pattern = Chem.MolFromSmarts('C1OC([H])C([H])C1O')  # Simplified sugar pattern

    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)O')

    # Check for sugar or phosphate groups
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if has_sugar:
        return False, "Contains sugar group; likely a nucleoside or nucleotide"
    if has_phosphate:
        return False, "Contains phosphate group; likely a nucleotide"

    # Check for nucleobase core structures
    for name, pattern in nucleobase_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains nucleobase core structure similar to {name}"

    # Check for modified nucleobase analogues using generic patterns
    pyrimidine_pattern = Chem.MolFromSmarts('c1cc[nH]c(=O)[nH]1')  # Generic pyrimidine core
    purine_pattern = Chem.MolFromSmarts('c1nc2c(n1)[nH]c(nc2)[N]')  # Generic purine core

    if mol.HasSubstructMatch(pyrimidine_pattern):
        return True, "Contains modified pyrimidine-like nucleobase core"
    if mol.HasSubstructMatch(purine_pattern):
        return True, "Contains modified purine-like nucleobase core"

    return False, "Does not contain nucleobase core structure"
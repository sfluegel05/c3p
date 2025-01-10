"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    
    has_carboxylic = mol.HasSubstructMatch(carboxylic_acid_pattern)
    has_carboxylate = mol.HasSubstructMatch(carboxylate_pattern)
    
    if not (has_carboxylic or has_carboxylate):
        return False, "No carboxylic acid group found"

    # Check for various types of amino groups
    amino_patterns = [
        Chem.MolFromSmarts("[NX3H2]"),  # Primary amine (-NH2)
        Chem.MolFromSmarts("[NX3H1]"),  # Secondary amine (-NH-)
        Chem.MolFromSmarts("[NX3H0]"),  # Tertiary amine (-N<)
        Chem.MolFromSmarts("[NX2]"),    # Imine (=NH)
        Chem.MolFromSmarts("[NH3+]"),   # Protonated amine
        Chem.MolFromSmarts("[NX3][H]")  # Explicit H notation
    ]
    
    has_amino = False
    for pattern in amino_patterns:
        if mol.HasSubstructMatch(pattern):
            has_amino = True
            break
            
    if not has_amino:
        return False, "No amino group found"

    # Additional checks to filter out some false positives
    
    # Check if the molecule is too small (less than 3 atoms)
    if mol.GetNumAtoms() < 3:
        return False, "Molecule too small to be an amino acid"
    
    # Check if it's just a simple amide (which has both COOH and NH2 but isn't an amino acid)
    amide_pattern = Chem.MolFromSmarts("[NX3H2]C(=O)[OH]")
    if mol.HasSubstructMatch(amide_pattern) and mol.GetNumAtoms() < 5:
        return False, "Simple amide structure"

    # If we've made it here, it's likely an amino acid
    reason = "Contains both carboxylic acid and amino groups"
    
    # Add detail about multiple amino groups if present
    amino_count = 0
    for pattern in amino_patterns:
        matches = mol.GetSubstructMatches(pattern)
        amino_count += len(matches)
    
    if amino_count > 1:
        reason += f" (found {amino_count} amino groups)"

    return True, reason
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon
    atoms and no glycosidic bonds to other sugar molecules.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms (C)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms for a monosaccharide"

    # Look for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found"

    # Look for carbonyl group (C=O)
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Check for glycosidic bonds (1->n connections)
    glycosidic_bond_pattern = Chem.MolFromSmarts("O-C-O")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        # Checking for potential larger sugar structure
        return False, "Potential oligosaccharide/polysaccharide structure detected"

    # Assess if molecule is likely a polyhydroxy aldehyde or ketone
    aldehyde_pattern = Chem.MolFromSmarts("O=CC(O)")
    ketone_pattern = Chem.MolFromSmarts("O=C[CH2]C(O)")
    if mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern):
        return True, "Structure matches polyhydroxy aldehyde/ketone pattern"
    
    return False, "Does not match monosaccharide structural requirements"
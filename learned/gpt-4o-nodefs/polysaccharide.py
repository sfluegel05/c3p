"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide has repeated units of monosaccharides linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a generic glycosidic linkage (e.g., -O-C-)
    glycosidic_pattern = Chem.MolFromSmarts("O-C")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
    
    # Look for multiple sugar ring units (more flexible pattern)
    cyclic_sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # Simplified monosaccharide pattern
    repeat_units = mol.GetSubstructMatches(cyclic_sugar_pattern)
    if len(repeat_units) < 2:
        return False, "Too few repeating monosaccharide units for polysaccharide"
    
    # Ensure the structure is large enough to be polysaccharide
    if len(mol.GetAtoms()) < 50:  # Arbitrarily chosen, can be adjusted
        return False, "Molecule too small for polysaccharide"
    
    # Polysaccharides should have extensive oxygen atoms due to many OH and glycosidic bonds
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 10:  # Arbitrarily chosen, can be adjusted based on known structures
        return False, "Too few oxygen atoms suggest an incomplete or incorrect polysaccharide"

    return True, "Contains multiple repeated monosaccharide units linked by glycosidic bonds"
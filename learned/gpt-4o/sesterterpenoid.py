"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid typically has a C25 skeleton derived from a sesterterpene,
    and may include rearrangements or modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms, sesterterpenoids are usually close to 25
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 22 or c_count > 35:
        return False, "Carbon count not typical for sesterterpenoids"

    # Improve isoprene unit detection angles 
    isoprene_pattern = Chem.MolFromSmarts("C(C)(C)C=C")  # consider other arrangements
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        # Less than 2 but matches complex patterns might be plausible
        five_membered_ring_pattern = Chem.MolFromSmarts("[R1]C1CCC1")  # cyclic terpenoid pattern
        if not mol.HasSubstructMatch(five_membered_ring_pattern):
            return False, "Insufficient or no isoprene unit or ring matches found"
    
    # Incorporate other functional group checks
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")  # simple ketone check
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(hydroxyl_matches) == 0 and len(ketone_matches) == 0:
        return False, "No typical terpenoid functional groups like hydroxyl or ketone found"
    
    # Recognize restructuring or common groups
    lactone_formation = Chem.MolFromSmarts("C1OC(=O)O1") # lactone ring common in some terpenoids
    if mol.HasSubstructMatch(lactone_formation):
        return True, "Lactone ring present, consistent with some sesterterpenoids"
    
    return True, "Molecule matches structural and functional criteria typical of sesterterpenoids"
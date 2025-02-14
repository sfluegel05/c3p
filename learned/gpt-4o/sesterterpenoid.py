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

    # Count carbon atoms, sesterterpenoids should approximate 25
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 35:
        return False, "Carbon count not typical for sesterterpenoids"
    
    # SMARTS patterns for identifying typical terpenoid units and structures
    isoprene_pattern = Chem.MolFromSmarts("C(C)(C)=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    five_membered_ring_pattern = Chem.MolFromSmarts("[R1]C1CCC1")
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)[CH2]O1")  # representative lactone rings
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")

    # Check for isoprene units, cyclic terpenoid patterns, or lactone features
    if (len(isoprene_matches) < 2 and 
        not mol.HasSubstructMatch(five_membered_ring_pattern) and
        not mol.HasSubstructMatch(lactone_pattern)):
        return False, "Insufficient structural motifs typical of sesterterpenoids"
    
    # Check for functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)[CX4]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if len(hydroxyl_matches) == 0 and len(ketone_matches) == 0 and len(ester_matches) == 0 and len(ether_matches) == 0:
        return False, "No relevant functional groups like hydroxyl, ketone, ester, or ether found"

    return True, "Molecule matches structural and functional criteria typical of sesterterpenoids"
"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:17915 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring compounds derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for prostanoic acid core structure with a more relaxed pattern
    prostanoic_pattern = Chem.MolFromSmarts("[C@@]12[C@@H]([C@@H]([C@@]1(CCCC(=O)[O,X1])[H])[H])[C@@H](C[C@@]2([H])[O,X1])[C,X2]=O")
    if not mol.HasSubstructMatch(prostanoic_pattern):
        return False, "Does not contain a recognizable prostanoic acid core structure"
    
    # Count rings - prostaglandins have 2 rings
    ring_info = mol.GetRingInfo()
    n_rings = len(set(x for y in ring_info.AtomRings() for x in y))
    if n_rings != 2:
        return False, f"Incorrect number of rings ({n_rings}), prostaglandins have 2"
    
    # Check for cis-double bonds in the cyclopentenone ring
    cis_ring_pattern = Chem.MolFromSmarts("[CX3]/C=C/[CX3]")
    cis_ring_matches = mol.GetSubstructMatches(cis_ring_pattern)
    if len(cis_ring_matches) < 1:
        return False, "Missing cis-double bond in cyclopentenone ring"
    
    # Check for carbonyl groups
    n_carbonyls = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0)
    if n_carbonyls < 2:
        return False, "Prostaglandins must have at least 2 carbonyl groups"
    
    # Check for hydroxy groups
    n_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if n_hydroxy < 1:
        return False, "Prostaglandins must have at least 1 hydroxy group"
    
    # All checks pass, classify as prostaglandin
    return True, "Contains a recognizable prostanoic acid core structure with characteristic ring systems, carbonyl groups, and hydroxy groups"
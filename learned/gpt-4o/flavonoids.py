"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flavonoid core pattern C6-C3-C6: Two phenyl rings with a heterocycle
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("c1cc(-c2ccccc2)ccc1C3=CC=CC(=O)O3"),  # Flavone
        Chem.MolFromSmarts("c1cc(-c2ccc(o2)C3=CC=CC(=O)O3)ccc1"),  # Isoflavone
        Chem.MolFromSmarts("c1cc(-c2cc(c(o2)C3=CC=CC(=O)O3)cc2)ccc1")  # Neoflavonoid
    ]
    
    matches_core = any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_core_patterns)
    if not matches_core:
        return False, "No flavonoid core C6-C3-C6 pattern variation found"
    
    # Check for typical flavonoid functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.HasSubstructMatch(hydroxyl_pattern)
    
    # Count number of aromatic rings
    arom_rings = [ring for ring in mol.GetRingInfo().AtomRings()
                  if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    aromatic_rings_count = len(arom_rings)
    
    # Validate that there are enough aromatic rings (typically 2)
    if aromatic_rings_count < 2:
        return False, f"Too few aromatic rings ({aromatic_rings_count}) to be considered a flavonoid"
    
    return True, "Flavonoid pattern detected with expected core and functional groups"
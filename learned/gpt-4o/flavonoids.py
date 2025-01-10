"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Expanded flavonoid core patterns including variations
    flavonoid_core_patterns = [
        Chem.MolFromSmarts("c1cc(-c2ccc(-c3ccco3)o2)ccc1"),  # General C6-C3-C6 pattern
        Chem.MolFromSmarts("c1cc(-c2ccccc2-c3ccc[oH]c3)occ1"),  # Isoflavonoids
        Chem.MolFromSmarts("c1ccc(cc1)-c2c[cH]cc(o2)c3cccc[oH]c3")  # Neoflavonoids
    ]

    matches_core = any(mol.HasSubstructMatch(pattern) for pattern in flavonoid_core_patterns)
    if not matches_core:
        return False, "No flavonoid core C6-C3-C6 pattern variation found"

    # Check for presence of carbonyl groups, typical in flavonoids
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Count number of aromatic rings (many flavonoids have high aromatic character)
    aromatic_rings_count = sum(1 for ring in mol.GetRingInfo().AtomRings()
                               if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if aromatic_rings_count < 2:
        return False, f"Too few aromatic rings ({aromatic_rings_count}) to be considered a flavonoid"

    return True, "Matches flavonoid pattern with expected core and functional groups resembling flavonoids"
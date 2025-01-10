"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid typically has a hydroxy-5beta-cholanic acid structure.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Bile acid steroid backbone pattern: cyclopenta[a]phenanthrene skeleton
    # with specific functional groups and stereochemistry.
    steroid_pattern = Chem.MolFromSmarts("C1C=C2[C@@H](CCC3)[C@H](C=C3)[C@@H]2C[C@@H]1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not match basic steroid backbone of bile acids"
    
    # Check for 5beta configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    # Here we simplify by checking that multiple chiral centers exist, though
    # proper 5beta stereochemistry is more complex to verify directly
    if len(chiral_centers) < 3:
        return False, f"Insufficient chiral centers for 5beta configuration, found {len(chiral_centers)}"
    
    # Look for common substitutions: hydroxyl groups and terminal carboxylic acid
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0)
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"
    
    carboxy_group = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxy_group):
        return False, "No terminal carboxylic acid group found"
    
    return True, "Matches bile acid structural criteria"
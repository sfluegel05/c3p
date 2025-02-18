"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:CHEBI_00000 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C15 skeleton of two benzene rings (A and B) connected by a three-carbon chain.
    Includes subclasses like flavones, isoflavones, chalcones, etc.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define core flavonoid skeleton patterns
    # Pattern 1: Two aromatic rings connected via a 3-carbon bridge (may include ketone for flavones)
    pattern1 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX3]-[CX3]=[OX1]")
    pattern2 = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX3H1]-[CX3H2]")  # For flavanones/chalcones
    
    # Check for matches to either core pattern
    match1 = mol.HasSubstructMatch(pattern1)
    match2 = mol.HasSubstructMatch(pattern2)
    
    if not (match1 or match2):
        return False, "No flavonoid core structure detected"
    
    # Check for presence of two aromatic rings (A and B)
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "Insufficient aromatic rings"
    
    # Basic carbon count check (C15 base)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Carbon count too low ({c_count})"
    
    return True, "Contains flavonoid core structure with two aromatic rings and three-carbon bridge"
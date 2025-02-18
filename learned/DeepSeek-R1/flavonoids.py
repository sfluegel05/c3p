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
    Flavonoids are characterized by a C15 skeleton of two aromatic rings (A and B) connected by a three-carbon bridge.
    Includes subclasses like flavones, isoflavones, chalcones, etc.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define core flavonoid patterns with more specificity
    # Pattern 1: Benzopyran core (flavones/flavonols) - two aromatic rings connected via oxygen-containing bridge
    benzopyran = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[#6]=[#6]-[#6](-[#8])=[#8]")
    
    # Pattern 2: Chalcone structure (open chain)
    chalcone = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[#6]=[#6]-[#6](=O)-[c]")
    
    # Pattern 3: Flavanone structure (saturated 3-carbon bridge)
    flavanone = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[#6]-[#6](-[#8])-[#6]-[c]")
    
    # Check for any core pattern match
    core_match = any(mol.HasSubstructMatch(p) for p in [benzopyran, chalcone, flavanone])
    
    if not core_match:
        return False, "No flavonoid core structure detected"
    
    # Verify presence of two aromatic rings connected via bridge
    # Get all aromatic rings
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "Insufficient aromatic rings"
    
    # Check if rings are connected through a 3-carbon bridge (including possible oxygen)
    # Find paths between aromatic rings that are 3-5 atoms long (allowing for O in bridge)
    ring_atoms = [set(ring) for ring in aromatic_rings]
    for i in range(len(ring_atoms)):
        for j in range(i+1, len(ring_atoms)):
            paths = Chem.GetShortestPaths(mol, ring_atoms[i].pop(), ring_atoms[j].pop(), maxLength=5)
            if any(3 <= len(path) <=5 for path in paths):
                return True, "Contains flavonoid core structure with two aromatic rings connected by three-carbon bridge"
    
    return False, "Aromatic rings not properly connected"
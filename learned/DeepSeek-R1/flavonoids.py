"""
Classifies: CHEBI:72544 flavonoids
"""
"""
Classifies: CHEBI:47916 flavonoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by two aromatic rings connected via a three-carbon bridge,
    which may be part of a heterocyclic ring. Includes subclasses like flavones, chalcones, etc.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define core flavonoid patterns with specific SMARTS
    patterns = {
        "flavone": Chem.MolFromSmarts("[c]1[c,c-,c=,c#][c,c-,c=,c#][c,c-,c=,c#][c,c-,c=,c#]1-[#6]@[#6](-[#8])=[#8]"),  # Benzopyrone core
        "chalcone": Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[#6]=[#6]-[#6](=[#8])-[c]"),  # Enone linker
        "isoflavone": Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[#6]1:[#6]:[#6]-[#8]-[#6]1=[#8]"),  # B-ring at position 3
        "flavanone": Chem.MolFromSmarts("[C&!a]1[C&!a][C&!a](=[O])[O][C&!a]1-[c]"),  # Saturated C-ring
        "aurone": Chem.MolFromSmarts("[O]=[C]1[C](=[O])[c][c][c][c]1-[O]-[c]")  # Benzofuranone
    }
    
    # Check for any core pattern match
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} core structure"
    
    # Fallback: Check aromatic rings and bridge connectivity
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() 
                      if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    
    if len(aromatic_rings) < 2:
        return False, "Insufficient aromatic rings"
    
    # Check all pairs of aromatic rings for bridge connectivity
    for i in range(len(aromatic_rings)):
        for j in range(i+1, len(aromatic_rings)):
            ring_a = set(aromatic_rings[i])
            ring_b = set(aromatic_rings[j])
            
            # Check shortest path between any two atoms from different rings
            for a1 in ring_a:
                for a2 in ring_b:
                    try:
                        path = rdmolops.GetShortestPath(mol, a1, a2)
                    except:
                        continue
                    if 3 <= len(path) <= 5:  # Allow 3-5 atom bridge (including possible ring bonds)
                        return True, "Two aromatic rings connected via three-carbon bridge"
    
    # Verify minimum size (C15 skeleton + substitutions)
    if rdMolDescriptors.CalcNumAtoms(mol) < 15:
        return False, "Too small for flavonoid skeleton"
    
    return False, "No characteristic flavonoid features detected"
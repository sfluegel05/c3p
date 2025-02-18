"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:15953 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide consists of two monosaccharide units joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    ring_info = mol.GetRingInfo()
    
    # Find all oxygen atoms not in any ring
    bridging_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        
        # Check if oxygen is in any ring
        in_ring = False
        for ring in ring_info.AtomRings():
            if atom.GetIdx() in ring:
                in_ring = True
                break
        if in_ring:
            continue
        
        # Check oxygen has exactly two neighbors (carbon or other)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            continue
        
        # Check if neighbors are in different connected components when oxygen is removed
        # Create a copy of the molecule without the oxygen
        modified_mol = Chem.RWMol(mol)
        modified_mol.RemoveAtom(atom.GetIdx())
        modified_mol = modified_mol.GetMol()
        
        # Get connected components
        parts = list(Chem.GetMolFrags(modified_mol, asMols=True))
        if len(parts) == 2:
            # Check if each part has characteristics of a monosaccharide
            # (at least 3 oxygen atoms as a heuristic)
            valid = True
            for part in parts:
                o_count = sum(1 for a in part.GetAtoms() if a.GetAtomicNum() == 8)
                if o_count < 3:
                    valid = False
                    break
            if valid:
                bridging_oxygens.append(atom.GetIdx())
    
    if len(bridging_oxygens) == 1:
        return True, "Two monosaccharide units connected by a glycosidic bond"
    else:
        return False, f"Found {len(bridging_oxygens)} valid glycosidic bonds, need exactly 1"
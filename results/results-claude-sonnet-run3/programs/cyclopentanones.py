from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclopentanones(smiles: str):
    """
    Determines if a molecule is a cyclopentanone (cyclopentane with at least one ketone group).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cyclopentanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for 5-membered rings
    rings = mol.GetRingInfo()
    five_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            five_rings.append(ring)
            
    if not five_rings:
        return False, "No 5-membered rings found"
        
    # For each 5-membered ring, check if it contains a ketone group
    for ring in five_rings:
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        
        # Check if ring is alicyclic (all sp3 carbons except ketone carbon)
        if any(atom.GetIsAromatic() for atom in ring_atoms):
            continue
            
        # Look for ketone pattern: C(=O) where C is in the ring
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if (neighbor.GetSymbol() == 'O' and 
                        neighbor.GetHybridization() == Chem.HybridizationType.SP2):
                        # Found ketone group in 5-membered ring
                        num_ketones = sum(1 for a in ring_atoms 
                                        for n in a.GetNeighbors()
                                        if n.GetSymbol() == 'O' and 
                                        n.GetHybridization() == Chem.HybridizationType.SP2)
                        return True, f"Cyclopentanone with {num_ketones} ketone group(s)"
                        
    return False, "No ketone group found in 5-membered alicyclic ring"
# Pr=1.0
# Recall=1.0
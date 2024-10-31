from rdkit import Chem
from rdkit.Chem import AllChem

def is_trithiane(smiles: str):
    """
    Determines if a molecule is a trithiane (ring with 3 sulfur and 3 carbon atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trithiane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()
    
    # Look for 5 or 6-membered rings
    valid_rings = []
    for ring in rings.AtomRings():
        if len(ring) in [5,6]:
            valid_rings.append(ring)
            
    if not valid_rings:
        return False, "No 5 or 6-membered rings found"

    # For each ring, check if it contains exactly 3 sulfur atoms
    for ring in valid_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        sulfur_count = sum(1 for atom in atoms if atom.GetSymbol() == 'S')
        carbon_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C')
        
        # Check atom counts
        if sulfur_count == 3 and carbon_count == len(ring)-3:
            # Check if ring is aliphatic (saturated)
            if all(not atom.GetIsAromatic() for atom in atoms):
                # Check if all carbons are sp3 hybridized
                if all(atom.GetHybridization() == Chem.HybridizationType.SP3 
                      for atom in atoms if atom.GetSymbol() == 'C'):
                    
                    # Get substituents
                    substituents = []
                    ring_atoms = set(ring)
                    for atom_idx in ring_atoms:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() not in ring_atoms:
                                substituents.append(neighbor.GetSymbol())
                    
                    if substituents:
                        return True, f"Trithiane with substituents: {', '.join(set(substituents))}"
                    else:
                        return True, "Unsubstituted trithiane"
                    
    return False, "No ring with exactly 3 sulfur atoms in aliphatic ring found"
# Pr=1.0
# Recall=1.0
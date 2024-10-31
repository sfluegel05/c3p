from rdkit import Chem
from rdkit.Chem import AllChem

def is_oxetanes(smiles: str):
    """
    Determines if a molecule is an oxetane or substituted oxetane.
    Oxetanes are 4-membered cyclic ethers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxetane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo()
    
    # Look for 4-membered rings
    four_membered_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 4:
            four_membered_rings.append(ring)
            
    if not four_membered_rings:
        return False, "No 4-membered rings found"

    # For each 4-membered ring, check if it contains exactly one oxygen
    for ring in four_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        oxygen_atoms = [atom for atom in atoms if atom.GetSymbol() == 'O']
        
        if len(oxygen_atoms) == 1:
            # Check if oxygen is sp3 hybridized
            if oxygen_atoms[0].GetHybridization() == Chem.HybridizationType.SP3:
                # Check if other atoms are carbons
                other_atoms = [atom for atom in atoms if atom.GetSymbol() != 'O']
                if all(atom.GetSymbol() == 'C' for atom in other_atoms):
                    # Count substituents
                    substituents = []
                    ring_atoms = set(ring)
                    
                    for atom_idx in ring_atoms:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIdx() not in ring_atoms:
                                substituents.append(neighbor.GetSymbol())
                    
                    if len(substituents) > 0:
                        return True, f"Substituted oxetane with substituents: {', '.join(set(substituents))}"
                    else:
                        return True, "Unsubstituted oxetane"
                    
    return False, "No oxetane ring found (4-membered ring with one sp3 oxygen and three carbons)"
# Pr=1.0
# Recall=0.875
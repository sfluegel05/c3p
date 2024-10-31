from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monobactam(smiles: str):
    """
    Determines if a molecule is a monobactam based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monobactam, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for presence of beta-lactam ring (4-membered ring with N and C=O)
    beta_lactam_pattern = Chem.MolFromSmarts('[N]1[C](=O)[C][C]1')
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
        
    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Get ring sizes
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    
    # Check for monocyclic structure (only one 4-membered ring)
    four_membered_rings = ring_sizes.count(4)
    if four_membered_rings != 1:
        return False, f"Found {four_membered_rings} 4-membered rings, expected 1"
        
    if num_rings > 1:
        # Allow for aromatic rings in substituents
        non_aromatic_rings = 0
        for ring in ring_info.AtomRings():
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if not all(atom.GetIsAromatic() for atom in atoms):
                non_aromatic_rings += 1
                
        if non_aromatic_rings > 1:
            return False, "Contains multiple non-aromatic rings"
    
    # Check for substituents on beta-lactam ring
    match = mol.GetSubstructMatch(beta_lactam_pattern)
    if not match:
        return False, "Could not map beta-lactam pattern"
        
    lactam_atoms = set(match)
    substituents = []
    
    for atom_idx in lactam_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in lactam_atoms:
                substituents.append(neighbor.GetSymbol())
                
    return True, f"Monobactam with beta-lactam ring and substituents: {', '.join(set(substituents))}"
# Pr=0.9473684210526315
# Recall=0.9
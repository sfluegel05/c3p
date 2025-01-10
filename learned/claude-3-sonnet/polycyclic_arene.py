"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a compound with multiple fused aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get ring information
    ri = mol.GetRingInfo()
    if ri.NumRings() < 2:
        return False, "Must contain at least 2 rings"
    
    # Check if molecule has aromatic atoms
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic rings found"
    
    # Get the aromatic rings
    aromatic_rings = []
    for ring in ri.AtomRings():
        # Check if ring is 6-membered and all atoms are aromatic
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(ring)
    
    if len(aromatic_rings) < 2:
        return False, "Must contain at least 2 aromatic 6-membered rings"
    
    # Check if rings are properly fused (share exactly 2 adjacent atoms)
    ring_sets = [set(ring) for ring in aromatic_rings]
    proper_fusion = False
    
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            shared_atoms = ring_sets[i].intersection(ring_sets[j])
            if len(shared_atoms) == 2:
                # Check if shared atoms are bonded
                shared_list = list(shared_atoms)
                if mol.GetBondBetweenAtoms(shared_list[0], shared_list[1]) is not None:
                    proper_fusion = True
                    break
        if proper_fusion:
            break
            
    if not proper_fusion:
        return False, "Aromatic rings must be properly fused (share an edge)"
    
    # Count carbon proportion - be more lenient to allow for substituents
    total_atoms = mol.GetNumAtoms()
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_atoms / total_atoms < 0.6:  # Reduced from 0.7
        return False, "Too few carbon atoms for a polycyclic arene"
    
    # Check if the core structure is primarily aromatic rings
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms / total_atoms < 0.5:  # Reduced from 0.6
        return False, "Core structure is not primarily aromatic"
    
    # Additional check for allowed substituents
    allowed_elements = {1, 6, 7, 8, 9, 17}  # H, C, N, O, F, Cl
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, "Contains disallowed elements"
            
    return True, "Contains multiple fused aromatic rings with primarily carbon atoms"
"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: polycyclic arene (polycyclic aromatic hydrocarbon)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDecomposition import GetSpanningTree

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
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(ring)
    
    if len(aromatic_rings) < 2:
        return False, "Must contain at least 2 aromatic rings"
    
    # Check if rings are fused
    # Convert ring atom lists to sets for intersection checks
    ring_sets = [set(ring) for ring in aromatic_rings]
    
    # Check if any two rings share atoms (are fused)
    fused_found = False
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            if ring_sets[i].intersection(ring_sets[j]):
                fused_found = True
                break
        if fused_found:
            break
            
    if not fused_found:
        return False, "Aromatic rings must be fused"
    
    # Count carbon proportion
    total_atoms = mol.GetNumAtoms()
    carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Most atoms should be carbon (allowing for some substituents)
    if carbon_atoms / total_atoms < 0.7:
        return False, "Too few carbon atoms for a polycyclic arene"
    
    # Check if the core structure is primarily aromatic rings
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms / total_atoms < 0.6:
        return False, "Core structure is not primarily aromatic"
        
    return True, "Contains multiple fused aromatic rings with primarily carbon atoms"
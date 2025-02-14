"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of fused aromatic rings.

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

    # Check for presence of aromatic rings
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic rings found"
    
    # Check if there are multiple aromatic rings
    aromatic_count = len(mol.GetSubstructMatches(aromatic_ring_pattern))
    if aromatic_count < 2:
        return False, "Less than 2 aromatic rings found, not a polycyclic arene"

    # Get ring information from RDKit
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    if not atom_rings:
        return False, "No rings detected"

    # Check for fused rings (at least two shared atoms)
    fused_rings = False
    num_rings = len(atom_rings)
    if num_rings > 1:
        for i in range(num_rings):
            for j in range(i + 1, num_rings):
                ring1 = set(atom_rings[i])
                ring2 = set(atom_rings[j])
                if len(ring1.intersection(ring2)) >= 2: # At least two shared atoms
                   fused_rings = True
                   break
            if fused_rings:
                break
    
    if not fused_rings:
       return False, "No fused aromatic rings detected"


    # Check for non-aromatic atoms (excluding hydrogens)
    non_aromatic_atoms = [atom for atom in mol.GetAtoms() if not atom.GetIsAromatic() and atom.GetAtomicNum() != 1]
    
    if len(non_aromatic_atoms) > mol.GetNumAtoms() * 0.3:
         return False, f"Too many non-aromatic atoms, {len(non_aromatic_atoms)} found."

    # Check that at least 50% of the atoms are aromatic.
    total_atoms = mol.GetNumAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())

    if total_atoms == 0:
        return False, "Molecule has no atoms"

    if aromatic_atoms/total_atoms < 0.5:
        return False, "Not enough aromatic atoms"

    return True, "Contains fused aromatic rings, classified as polycyclic arene"
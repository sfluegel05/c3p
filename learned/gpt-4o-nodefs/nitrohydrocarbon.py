"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon contains nitro groups (-NO2) attached to a primarily hydrocarbon framework.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)

    if len(nitro_matches) == 0:
        return False, "No nitro groups found"

    # Check the carbon framework
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    total_atom_count = mol.GetNumAtoms()
    o_n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8])

    if c_count / total_atom_count < 0.5:
        return False, "Insufficient carbon atoms in framework"

    # Check that the molecule doesn't have excessive heteroatoms excluding nitro
    non_c_h_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 7, 8]]
    if len(non_c_h_atoms) > 0:
        return False, "Contains heteroatoms not typical in hydrocarbon framework"

    # Check for the dominance of carbon in ring systems
    ri = mol.GetRingInfo()
    if ri.NumRings() > 0:
        for bond_ring in ri.BondRings():
            atoms_in_ring = {mol.GetBondWithIdx(bidx).GetBeginAtomIdx() for bidx in bond_ring}
            carbons_in_ring = sum(1 for idx in atoms_in_ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            
            # Ensure rings are primarily carbon-based
            if carbons_in_ring / len(atoms_in_ring) < 0.6:
                return False, "Ring structures present, not predominantly hydrocarbon"

    return True, "Contains nitro groups attached to a primarily hydrocarbon framework"
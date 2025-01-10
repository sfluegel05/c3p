"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:8338 porphyrins
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a fundamental skeleton of four pyrrole nuclei united through the alpha-positions by four methine groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure aromaticity and ring information is computed
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.MolSanitizeException as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Search for 16-membered rings
    sixteen_membered_rings = [ring for ring in atom_rings if len(ring) == 16]
    if not sixteen_membered_rings:
        return False, "No 16-membered ring found"

    # For each 16-membered ring, check for porphyrin characteristics
    for ring in sixteen_membered_rings:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Count nitrogen atoms in the ring
        nitrogen_atoms = [atom for atom in atoms_in_ring if atom.GetAtomicNum() == 7]
        if len(nitrogen_atoms) != 4:
            continue  # Not a porphyrin if not exactly 4 nitrogen atoms
        
        # Check that nitrogens are pyrrolic ([nH])
        pyrrolic_nitrogens = [
            atom for atom in nitrogen_atoms
            if atom.GetIsAromatic() and atom.GetTotalNumHs() == 1
        ]
        if len(pyrrolic_nitrogens) != 4:
            continue  # Not all nitrogens are pyrrolic

        # Check that the ring is aromatic
        aromatic_atoms = [atom for atom in atoms_in_ring if atom.GetIsAromatic()]
        if len(aromatic_atoms) != len(atoms_in_ring):
            continue  # Ring is not fully aromatic

        # Passed all checks; porphyrin ring system found
        return True, "Porphyrin ring system detected"

    # If no rings passed the checks
    return False, "No porphyrin ring system detected"
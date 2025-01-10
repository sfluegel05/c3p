"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam features a four-membered ring that includes the amide nitrogen
    and the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the beta-lactam core
    # '[R2]1C(=O)NC1' represents a four-membered ring (R2) with C(=O) and N
    beta_lactam_pattern = Chem.MolFromSmarts("[R2]1C(=O)NC1")
    if mol.HasSubstructMatch(beta_lactam_pattern):
        # Verify that the matched pattern corresponds to a ring
        ring_info = mol.GetRingInfo()
        for ring_atoms in ring_info.AtomRings():
            if len(ring_atoms) == 4:
                # A valid beta-lactam must have C(=O) and N in the ring atoms
                matching_atom_indices = mol.GetSubstructMatch(beta_lactam_pattern)
                if set(matching_atom_indices).issubset(set(ring_atoms)):
                    return True, "Contains a four-membered beta-lactam ring"
    
    return False, "Does not contain a four-membered beta-lactam ring"
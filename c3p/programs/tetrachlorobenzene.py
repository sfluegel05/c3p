"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying four chloro groups, where the benzene ring (c1ccccc1) is substituted with exactly 4 chlorine atoms attached.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    The algorithm looks for aromatic benzene rings (6 carbons, aromatic) and checks whether any one of these rings
    has exactly four chlorine atoms (atomic number 17) directly attached as substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a benzene ring substituted with exactly 4 Cl atoms, False otherwise.
        str: The reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure aromaticity is perceived
    Chem.SanitizeMol(mol)
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"

    # Loop over each ring in the molecule
    for ring in ring_info.AtomRings():
        # Check if the ring has 6 atoms (benzene candidate)
        if len(ring) != 6:
            continue
        
        # Verify all ring atoms are carbon and are aromatic 
        is_benzene_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene_ring = False
                break
        if not is_benzene_ring:
            continue

        # Check the number of chlorine substituents attached to this benzene ring.
        # For each ring atom, check neighbors that are NOT part of the ring and if they are Cl.
        cl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # If the neighbor is not in the ring, consider it a substituent
                if neighbor.GetIdx() not in ring:
                    if neighbor.GetAtomicNum() == 17:  # Chlorine atomic number
                        cl_count += 1
        if cl_count == 4:
            return True, "Found benzene ring with exactly 4 chlorine substituents"
    
    # If no benzene ring with exactly 4 Cl substituents is found, return False.
    return False, "No benzene ring with exactly 4 chlorine substituents found"
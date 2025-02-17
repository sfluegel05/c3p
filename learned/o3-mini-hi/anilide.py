"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide - Any aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem

def is_in_benzene_ring(atom, mol):
    """
    Check if the given atom is part of a benzene ring.
    A benzene ring is defined as a six-membered ring where every atom is a carbon and aromatic.
    """
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    # Iterate over each ring (as a tuple of atom indices)
    for ring in ring_info.AtomRings():
        if len(ring) == 6 and atom.GetIdx() in ring:
            # Check that every atom in the ring is a carbon and is aromatic.
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                return True
    return False

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as an aromatic amide where the nitrogen of the amide bond
    is derived from aniline. In other words, the amide nitrogen is directly attached to a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an anilide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find amide substructure matches.
    # We use a simple SMARTS for an amide: nitrogen (N) bonded to a carbon (C) that is double-bonded to oxygen.
    amide_smarts = "[NX3][CX3](=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    matches = mol.GetSubstructMatches(amide_query)
    
    if not matches:
        return False, "No amide functionality found"
    
    # For each match, check if the amide nitrogen has a substituent that is directly part of a benzene ring.
    for match in matches:
        n_idx = match[0]  # amide nitrogen index
        c_idx = match[1]  # carbonyl carbon index, ignore it for checking substituents
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check each neighbor of the amide nitrogen (except the carbonyl carbon)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == c_idx:
                continue  # skip the carbonyl carbon
            
            # The neighbor should be a carbon atom.
            if neighbor.GetAtomicNum() != 6:
                continue
            
            # Additionally, check that this carbon is part of a benzene ring
            if is_in_benzene_ring(neighbor, mol):
                return True, "Found amide group with the nitrogen directly attached to a benzene ring (anilide)"
    
    return False, "Amide group(s) found but none have the nitrogen attached directly to a benzene ring"

# Example usage (you can test with one of the provided SMILES strings):
if __name__ == "__main__":
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"  
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)
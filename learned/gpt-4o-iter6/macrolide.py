"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is characterized by a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has initialized ring information properly
    ring_info = mol.GetRingInfo()
    if not ring_info:
        return False, "Molecule ring information could not be extracted"

    # Define a potentially more flexible SMARTS pattern for macrolactone
    # We need a cyclic ester structure; allowing for possible variations in attached atoms
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)[C@@H]1")  # Extending this pattern can be necessary

    # Ensure ring info is correctly initialized
    if not ring_info.NumAtomRings():
        return False, "No rings identified in the molecule"

    # Iterate over all rings
    for ring_atoms in ring_info.AtomRings():
        if len(ring_atoms) >= 12:  # Check if ring size is 12 or more
            try:
                submol = Chem.PathToSubmol(mol, ring_atoms)
                if submol.HasSubstructMatch(lactone_pattern):  # Looking for lactone pattern
                    return True, "Contains a macrocyclic lactone with 12 or more members"
            except Exception as e:
                return False, f"Error in substructure matching: {str(e)}"
    
    return False, "Does not contain a macrocyclic lactone with 12 or more members"

# Example usage (uncomment to test):
# result, reason = is_macrolide("O=C1O[C@H](C(=O)NC[C@@H](O)C[C@@H](O)[C@@H](O)CCC(=CC=C[C@H](CCCC(CCCC[C@H](C[C@H]1C)C)=O)CC)COC)CCC")
# print(result, reason)
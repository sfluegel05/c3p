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
    
    # Define the SMARTS pattern for a lactone
    # Cycloalkane with a ester functionality needed
    lactone_pattern = Chem.MolFromSmarts("C1C(=O)O[C;R1]1")
    
    # Use RDKit's built-in methods to extract rings
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if len(ring_atoms) >= 12:  # Ensure ring size is 12 or more
            # Extract substructure corresponding to the ring
            submol = Chem.PathToSubmol(mol, ring_atoms)
            # Check for lactone pattern within this macroring
            if submol.HasSubstructMatch(lactone_pattern):
                return True, "Contains a macrocyclic lactone with 12 or more members"
    
    return False, "Does not contain a macrocyclic lactone with 12 or more members"

# Example usage:
# result, reason = is_macrolide("O=C1O[C@H](C(=O)NC[C@@H](O)C[C@@H](O)[C@@H](O)CCC(=CC=C[C@H](CCCC(CCCC[C@H](C[C@H]1C)C)=O)CC)COC)CCC")
# print(result, reason)
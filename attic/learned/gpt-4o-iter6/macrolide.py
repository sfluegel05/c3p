"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for a macrocyclic lactone (Cyclic ester with ring size of 12 or more)
    # SMARTS for esters: [C](=O)O
    ester_pattern = Chem.MolFromSmarts("[C](=O)O")
    # Get rings, filter for those containing an ester
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if len(ring_atoms) >= 12:
            submol = Chem.PathToSubmol(mol, ring_atoms)
            if submol.HasSubstructMatch(ester_pattern):
                return True, "Contains a macrocyclic lactone with 12 or more atoms"

    return False, "Does not contain a macrocyclic lactone with 12 or more atoms"


# Example usage:
# result, reason = is_macrolide("SMILES_STRING_HERE")
# print(result, reason)
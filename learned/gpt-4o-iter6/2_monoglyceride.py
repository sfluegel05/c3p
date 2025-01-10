"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is defined by the ester linkage at the second carbon of a glycerol backbone,
    leaving the first and third carbons with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a 2-monoglyceride
    # Glycerol backbone with the ester linkage at C2 and hydroxyl groups at the terminal positions:
    # [OH]-C([OH])(CO)-C(=O)
    glycerol_2_mono_pattern = Chem.MolFromSmarts("O(C(CO)CO)C(=O)")
    
    # Search for matches in the molecule against the pattern
    matches = mol.GetSubstructMatches(glycerol_2_mono_pattern)
    
    if not matches:
        return False, "No 2-monoglyceride pattern found with esterification at C2"

    # Iterate through the matches, confirming the ester linkage and hydroxyl placements
    for match in matches:
        c1, c2, ester_o = match[:3]

        # Verify C1 and C3 neighboring patterns to ensure hydroxyl presence
        atom_c1 = mol.GetAtomWithIdx(c1)
        atom_c2 = mol.GetAtomWithIdx(c2)

        # Ensure C1 and C3 have appropriate hydroxyl neighbors
        if not any(nbr.GetAtomicNum() == 8 for nbr in atom_c1.GetNeighbors()):
            return False, "C1 does not have a hydroxyl group"
        if not any(nbr.GetAtomicNum() == 8 for nbr in atom_c2.GetNeighbors()):
            return False, "C3 does not have a hydroxyl group"

        # Validate ester linkage correctly at C2
        ester_ox = mol.GetAtomWithIdx(ester_o)
        if not any(nbr.GetIdx() == c2 for nbr in ester_ox.GetNeighbors()):
            return False, "Ester linkage not correctly positioned at C2"

    return True, "Contains the structure of 2-monoglyceride with ester at C2"
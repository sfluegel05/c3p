"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more flexible SMARTS pattern for 1-benzopyran
    benzopyran_patterns = [
        Chem.MolFromSmarts("c1cc2occcc2c1"), # General 1-benzopyran
        Chem.MolFromSmarts("c1cc2[cH1]occc2c1"), # Allow for aromatic and hydrogen-bonding
        Chem.MolFromSmarts("C1=COc2ccccc2C1"), # Consider benzopyran with variation
    ]
    
    # Check for presence of 1-benzopyran structure
    for pattern in benzopyran_patterns:
        if mol.HasSubstructMatch(pattern):
            # Identify atoms that match the 1-benzopyran pattern
            benzopyran_matches = mol.GetSubstructMatches(pattern)

            # For each match, check if there's an aryl ring (aromatic six-membered ring) at position 4
            for match in benzopyran_matches:
                possible_aryl_positions = [4, 5]  # Possible indexing changes in different substructure mappings
                for pos_idx in possible_aryl_positions:
                    aryl_position_atom = mol.GetAtomWithIdx(match[pos_idx])
                    
                    # Scan the neighbors for an aryl group
                    for neighbor in aryl_position_atom.GetNeighbors():
                        # Check if this neighbor forms part of an aromatic six-membered ring
                        if neighbor.GetIsAromatic() and mol.GetRingInfo().IsAtomInRingOfSize(neighbor.GetIdx(), 6):
                            return True, "Aryl substituent found at or near position 4"
    
    return False, "No aryl substitution at or near position 4 in the 1-benzopyran"
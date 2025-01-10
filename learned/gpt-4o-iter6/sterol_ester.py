"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is obtained by formal condensation of 
    the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general sterol pattern (fused four-ring steroid structure)
    sterol_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC(O4)C")  # simplified four-ring structure

    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No typical sterol backbone found"

    # Define ester linkage pattern
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check if the ester linkage is part of the molecule
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    if not ester_matches:
        return False, "No ester linkage found"

    # Assume the linkage is to the 3-hydroxy group of sterol
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    # Check connectivity: ester linkage to the sterol hydroxyl
    for ester_match in ester_matches:
        ester_oxygen_idx = ester_match[2]  # Index of the oxygen in the ester
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        ester_neighbors = ester_oxygen_atom.GetNeighbors()
        
        # Check if any neighbor of the ester oxygen is the hydroxyl group (on correct carbon)
        for neighbor in ester_neighbors:
            if neighbor.HasSubstructMatch(hydroxyl_pattern):
                sterol_linked = False
                # Ensure that there's correct bonding for the linkage
                for sterol_atom_idx in sterol_pattern:
                    if mol.GetAtomWithIdx(sterol_atom_idx).GetSymbol() == 'C':
                        if est_oxygen_idx in [n.GetIdx() for n in mol.GetAtomWithIdx(sterol_atom_idx).GetNeighbors()]:
                            sterol_linked = True
                            break
                if sterol_linked:
                    return True, "Contains sterol backbone with ester linkage at the 3-hydroxy group"
    
    return False, "No ester linkage connected to the 3-hydroxy group of sterol found"
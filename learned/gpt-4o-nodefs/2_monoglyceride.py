"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride consists of a glycerol backbone with a single fatty acid attached
    at the 2 position (middle) via an ester bond, and hydroxyl groups on the 1 and 3 positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for a glycerol 2-mono ester: [C(OH)(CO)]
    # Glycerol backbone (with ester bond on 2 position): 1,3-dihydroxypropan-2-yl ester
    glycerol_pattern = Chem.MolFromSmarts("C(CO)(CO)OC(=O)[C;H0]")
    
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 2-monoglyceride structure found"
    
    # Check if this match is consistent for a 2-position attachment
    matches = mol.GetSubstructMatches(glycerol_pattern)
    for match in matches:
        attachment_index = match[2]  # 2-carbon position in glycerol
        ester_group = mol.GetAtomWithIdx(attachment_index)
        
        # Ensure the ester bond is at the correct position
        connected_atoms = [neighbor.GetSymbol() for neighbor in ester_group.GetNeighbors()]
        
        if 'O' not in connected_atoms:  # Bonding to oxygen in ester
            continue
            
        # Check other carbons in the glycerol backbone are hydroxylated
        hydroxyl_atoms = []
        for i in match:
            atom = mol.GetAtomWithIdx(i)
            if atom.GetSymbol() == 'C':
                neighbors = [neigh.GetSymbol() for neigh in atom.GetNeighbors()]
                if 'O' in neighbors:
                    hydroxyl_atoms.append('O')
        
        if len(hydroxyl_atoms) == 2:
            return True, "Contains glycerol backbone with fatty acid esterified at 2 position"
    
    return False, "Structure does not match 2-monoglyceride pattern"
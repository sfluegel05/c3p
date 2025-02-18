"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:39143 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has a nitrogen atom with +1 charge bonded to exactly two carbons (protonated secondary amine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        # Check for positively charged nitrogen
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Count carbon neighbors
            carbon_count = sum(1 for neighbor in atom.GetNeighbors() 
                              if neighbor.GetAtomicNum() == 6)
            
            # Secondary ammonium has exactly two carbon neighbors
            if carbon_count == 2:
                return True, "Positively charged nitrogen with two carbon substituents"
    
    return False, "No secondary ammonium group (N+ with two carbons) found"
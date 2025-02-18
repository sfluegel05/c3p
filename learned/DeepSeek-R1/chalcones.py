"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    Chalcones are 1,3-diphenylpropenones (Ar-CH=CH-CO-Ar) and derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define the core chalcone pattern: two aromatic rings connected by propenone (C=C-C=O)
    chalcone_pattern = Chem.MolFromSmarts("[c]C(=O)C=C[c]")
    matches = mol.GetSubstructMatches(chalcone_pattern)
    
    if not matches:
        return False, "No chalcone core (Ar-C(=O)-C=C-Ar) found"
    
    # Verify the carbonyl is a ketone (not part of ester/acid/amide)
    for match in matches:
        carbonyl_idx = match[1]  # Index of the carbonyl carbon in the SMARTS pattern
        atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Check neighbors of carbonyl carbon
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:  # Oxygen without H (ester/acid)
                return False, "Carbonyl group is part of an ester or acid"
            if neighbor.GetAtomicNum() == 7:  # Nitrogen (amide)
                return False, "Carbonyl group is part of an amide"
    
    return True, "Contains chalcone core (Ar-C(=O)-C=C-Ar) with ketone group"
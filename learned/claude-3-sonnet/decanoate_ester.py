"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:36392 decanoate ester
A fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) 
with the hydroxy group of an alcohol or phenol.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for decanoate (decanoyl) group pattern
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)[O;X2]")
    decanoate_matches = mol.GetSubstructMatches(decanoate_pattern)
    
    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check if decanoate group is present and attached via ester bond
    if len(decanoate_matches) > 0 and len(ester_matches) > 0:
        # Check if decanoate group is attached to ester oxygen
        for ester_match in ester_matches:
            ester_oxygen = ester_match[0]
            for decanoate_match in decanoate_matches:
                if ester_oxygen in decanoate_match:
                    return True, "Contains a decanoate (decanoyl) group attached via an ester bond"
    
    # If no decanoate ester group found
    return False, "Does not contain a decanoate (decanoyl) ester group"
"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
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
    
    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Pattern for decanoyl group (10-carbon chain with carbonyl)
    # Note: The pattern looks for C(=O)-CCCCCCCCC
    decanoyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    decanoyl_matches = mol.GetSubstructMatches(decanoyl_pattern)
    
    if not decanoyl_matches:
        return False, "No decanoyl group found"
    
    # Check if any decanoyl group is part of an ester
    found_decanoate_ester = False
    for ester_match in ester_matches:
        for decanoyl_match in decanoyl_matches:
            # Check if the carbonyl carbon of the ester matches the carbonyl carbon of the decanoyl
            if ester_match[1] == decanoyl_match[0]:
                found_decanoate_ester = True
                break
        if found_decanoate_ester:
            break
            
    if not found_decanoate_ester:
        return False, "Decanoyl group not connected via ester linkage"
    
    return True, "Contains decanoyl group connected via ester linkage"
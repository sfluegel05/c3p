"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for butyrate ester pattern: -O-C(=O)-CCC-C
    # This pattern allows for substitutions on the carbons
    butyrate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4][CX4][CX4]")
    
    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
        return False, "No butyrate chain found"
        
    # For each potential match, verify basic connectivity
    for match in matches:
        ester_o = mol.GetAtomWithIdx(match[0])
        carbonyl_c = mol.GetAtomWithIdx(match[1])
        c1 = mol.GetAtomWithIdx(match[2])
        c2 = mol.GetAtomWithIdx(match[3])
        c3 = mol.GetAtomWithIdx(match[4])
        
        # Verify we have carbons in the chain
        if all(atom.GetAtomicNum() == 6 for atom in [carbonyl_c, c1, c2, c3]):
            # Check that the ester oxygen is connected to another carbon
            # (to confirm it's an ester and not a carboxylic acid)
            for neighbor in ester_o.GetNeighbors():
                if neighbor.GetIdx() != carbonyl_c.GetIdx() and neighbor.GetAtomicNum() == 6:
                    return True, "Contains butyrate ester group (R-O-C(=O)-CH2-CH2-CH3)"
    
    return False, "Found similar structure but not a butyrate ester"
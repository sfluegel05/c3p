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

    # Look for basic ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Look for butyrate ester pattern with more flexible matching
    # Allow for potential substitutions on the carbon chain
    butyrate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[$([CX4H2]),$([CX4H1]),$([CX4H0])][$([CX4H2]),$([CX4H1]),$([CX4H0])][CX4]")
    matches = mol.GetSubstructMatches(butyrate_pattern)
    
    if not matches:
        return False, "No butyrate ester group found"
    
    # For each potential match, verify it's a butyrate ester
    for match in matches:
        if len(match) < 5:  # Safety check
            continue
            
        # Get the atoms involved in the match
        ester_o = mol.GetAtomWithIdx(match[0])
        carbonyl_c = mol.GetAtomWithIdx(match[1])
        c1 = mol.GetAtomWithIdx(match[2])
        c2 = mol.GetAtomWithIdx(match[3])
        c3 = mol.GetAtomWithIdx(match[4])
        
        # Verify the basic connectivity
        if not all([atom.GetDegree() <= 4 for atom in [c1, c2, c3]]):
            continue
            
        # Check that we have a proper carbon chain
        # Allow for substitutions but verify basic carbon backbone
        if (c1.GetAtomicNum() == 6 and 
            c2.GetAtomicNum() == 6 and 
            c3.GetAtomicNum() == 6):
            
            # Check that the chain is properly terminated
            # The last carbon should have at least one hydrogen
            if c3.GetTotalNumHs() >= 1:
                return True, "Contains butyrate ester group (R-O-C(=O)-CH2-CH2-CH3 or substituted variant)"
    
    return False, "Found similar structure but not a true butyrate ester"
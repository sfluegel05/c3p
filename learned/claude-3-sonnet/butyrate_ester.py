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

    # Look for ester group connected to a 3-carbon chain ending in methyl
    # [OX2] - ester oxygen
    # [CX3](=[OX1]) - carbonyl carbon with double-bonded oxygen
    # [C] - any carbon (can be substituted)
    # [C] - any carbon (can be substituted)
    # [CH3,CH2] - terminal carbon (methyl or methylene)
    butyrate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[C][C][CH3,CH2]")
    
    matches = mol.GetSubstructMatches(butyrate_pattern)
    if not matches:
        return False, "No butyrate ester group found"
        
    # For each potential match, verify it's a proper ester
    for match in matches:
        ester_o = mol.GetAtomWithIdx(match[0])
        carbonyl_c = mol.GetAtomWithIdx(match[1])
        
        # Verify the ester oxygen is connected to a carbon (making it an ester)
        for neighbor in ester_o.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_c.GetIdx() and neighbor.GetAtomicNum() == 6:
                # Verify the carbons in the chain are connected properly
                c1 = mol.GetAtomWithIdx(match[2])
                c2 = mol.GetAtomWithIdx(match[3])
                c3 = mol.GetAtomWithIdx(match[4])
                
                # Check if atoms form a continuous chain
                if (c1 in carbonyl_c.GetNeighbors() and 
                    c2 in c1.GetNeighbors() and 
                    c3 in c2.GetNeighbors()):
                    return True, "Contains butyrate ester group (R-O-C(=O)-CH2-CH2-CH3)"
    
    return False, "Structure contains ester but not a butyrate ester"
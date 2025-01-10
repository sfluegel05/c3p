"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:33566 aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are compounds with the general structure R-CH=N-OH, where R is any carbon group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Various SMARTS patterns for different oxime representations
    patterns = [
        # Standard aldoxime pattern: carbon must have one H and one single bond to carbon
        "[C!R][CH1](=[NX2]-[OX2H1])",  # R-CH=N-OH
        "[C!R][CH1](=[NX2]-[OX1-])",   # R-CH=N-O- (deprotonated)
        # E/Z configurations
        "[C!R][CH1](/=[NX2]/[OX2H1])",  # E isomer
        "[C!R][CH1](\\=[NX2]\\[OX2H1])", # Z isomer
        "[C!R][CH1](/=[NX2]\\[OX2H1])",  # Alternative E/Z
        "[C!R][CH1](\\=[NX2]/[OX2H1])"   # Alternative E/Z
    ]
    
    found_match = False
    for pattern in patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(smarts):
            matches = mol.GetSubstructMatches(smarts)
            for match in matches:
                # Get the carbon atom of the C=N bond
                carbon_idx = match[1]  # Index 1 is the CH1 atom
                carbon = mol.GetAtomWithIdx(carbon_idx)
                
                # Verify carbon has exactly one hydrogen
                if carbon.GetTotalNumHs() != 1:
                    continue
                    
                # Verify carbon is not aromatic
                if carbon.GetIsAromatic():
                    continue
                    
                # Verify carbon is not part of a ring
                if carbon.IsInRing():
                    continue
                    
                # Count number of single bonds to carbon atoms
                carbon_neighbors = [n for n in carbon.GetNeighbors() 
                                 if n.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(carbon_idx, n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
                
                if len(carbon_neighbors) == 1:
                    found_match = True
                    break
                    
            if found_match:
                break
    
    if not found_match:
        return False, "No aldoxime group (R-CH=N-OH) found"
        
    # Additional check for oxime esters
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)-[OX2]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
        for match in matches:
            oxime_carbon = match[1]
            for carbonyl_match in carbonyl_matches:
                if oxime_carbon == carbonyl_match[0]:
                    return False, "Contains oxime ester group instead of aldoxime"
    
    return True, "Contains aldoxime group (R-CH=N-OH)"
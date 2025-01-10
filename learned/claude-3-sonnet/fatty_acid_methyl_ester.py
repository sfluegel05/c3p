"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
A fatty acid ester that is the carboxylic ester obtained by the formal condensation 
of a fatty acid with methanol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for methyl ester group pattern
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if not methyl_ester_matches:
        return False, "No methyl ester group found"
    
    # Count number of methyl ester groups
    num_methyl_esters = len(methyl_ester_matches)
    
    # Special case: Allow dimethyl esters (like dimethyl sebacate)
    if num_methyl_esters > 2:
        return False, f"Too many methyl ester groups ({num_methyl_esters})"
    
    # Get atoms in methyl ester group
    ester_atoms = set()
    for match in methyl_ester_matches:
        ester_atoms.update(match)
    
    # Check for carbon chain attached to ester group
    carbon_chain = False
    for match in methyl_ester_matches:
        carbonyl_carbon = match[0]  # First atom in pattern is carbonyl carbon
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetIdx() not in ester_atoms and neighbor.GetAtomicNum() == 6:
                # Found carbon attached to ester that's not part of ester group
                carbon_chain = True
                break
    
    if not carbon_chain:
        return False, "No carbon chain attached to ester group"
    
    # Count carbons (excluding methyl ester carbons)
    non_ester_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in ester_atoms:
            non_ester_carbons += 1
            
    if non_ester_carbons < 2:
        return False, "Carbon chain too short for fatty acid"
        
    # If we have exactly two methyl esters, check if it's a valid diester
    if num_methyl_esters == 2:
        # Check for linear chain between esters
        chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4].[CX3](=[OX1])[OX2][CH3X4]")
        if not mol.HasSubstructMatch(chain_pattern):
            return False, "Invalid diester structure"
        return True, "Valid dimethyl ester of dicarboxylic acid"
    
    return True, "Contains methyl ester group with appropriate carbon chain"
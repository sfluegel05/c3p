"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: CHEBI:25495 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where one hydroxyl group is esterified with phosphoric acid and the other two with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern with consideration of stereochemistry and ester sites
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    
    # Define phosphoric acid ester pattern
    phosphoester_pattern = Chem.MolFromSmarts("COP(O)(O)=O")
    
    # Phosphatidic acid should include both the glycerol backbone and the phosphoric acid ester group
    if not mol.HasSubstructMatch(glycerol_pattern) or not mol.HasSubstructMatch(phosphoester_pattern):
        return False, "Lacks necessary glycerol backbone or phosphoric ester group"
    
    # Check for two esterified fatty acid groups with a defined length (at least 12 carbons typically)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    fatty_acid_count = 0
    for match in ester_matches:
        chain_length = len(match) - 2  # Compensating for CO group
        atom_count = 0
        for atom in match:
            if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6: # Counting carbons
                atom_count += 1
        if atom_count >= 12:
            fatty_acid_count += 1

    if fatty_acid_count < 2:
        return False, "Less than 2 long-chain fatty acid ester groups found"

    return True, "Contains glycerol backbone with 2 long-chain fatty acid esters and a phosphoric acid ester"
"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester should have lauric acid (12-carbon saturated chain) as the carboxylic acid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a dodecanoate ester
    lauric_acid_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O")
    if lauric_acid_ester_pattern is None:
        return None, None  # If the pattern cannot be compiled

    # Check for match with the pattern
    if mol.HasSubstructMatch(lauric_acid_ester_pattern):
        # Additionally, count carbons in the match to ensure a 12-carbon chain is present in the ester
        matches = mol.GetSubstructMatches(lauric_acid_ester_pattern)
        for match in matches:
            # Extract atoms involved in the ester linkage
            ester_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            # Count carbons in the alkyl chain
            carbon_count = sum(1 for atom in ester_atoms if atom.GetAtomicNum() == 6)
            if carbon_count == 12:
                return True, "Contains lauric acid (12-carbon) ester component"
        
    return False, "Does not contain lauric acid ester component"
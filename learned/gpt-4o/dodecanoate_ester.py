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

    # SMARTS pattern for a dodecanoate ester (12-carbon chain attached to ester group)
    # This pattern specifically checks for a C(=O)O linkage with a long alkyl chain
    lauric_acid_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O")  # 12 carbon chain before ester_group

    # Attempt to match the pattern
    if mol.HasSubstructMatch(lauric_acid_ester_pattern):
        # Verify the context of the matching substructure
        matches = mol.GetSubstructMatches(lauric_acid_ester_pattern)
        for match in matches:
            # Verify the chain length (encapsulating additional context checks if needed)
            carbon_atoms = [mol.GetAtomWithIdx(idx) for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
            if len(carbon_atoms) == 12:
                return True, "Contains lauric acid (12-carbon) ester component"
        
    return False, "Does not contain lauric acid ester component"
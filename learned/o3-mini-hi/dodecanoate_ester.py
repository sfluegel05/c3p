"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoate ester)
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester.
    A dodecanoate ester is defined as any fatty acid ester in which the carboxylic acid component is lauric acid 
    (i.e. the acyl group corresponds to CH3(CH2)10C(=O)-).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains the dodecanoate (laurate) ester moiety, False otherwise.
        str: Explanation for classification.
    """
    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define SMARTS for the laurate (dodecanoate) ester moiety.
    # Here "CCCCCCCCCCCC" represents a chain of 12 carbons (lauric acid in its acid form) 
    # followed by (C(=O)O) which represents the carbonyl-ester functionality.
    dodecanoate_smarts = "CCCCCCCCCCCC(=O)O"
    query = Chem.MolFromSmarts(dodecanoate_smarts)
    if query is None:
        return False, "Failed to create SMARTS pattern for dodecanoate ester."
    
    # Check if the molecule has a substructure match for the dodecanoate ester.
    if not mol.HasSubstructMatch(query):
        return False, "Dodecanoate ester substructure not found."
    
    return True, "Molecule contains a dodecanoate ester (laurate ester) group."

# Example usage (you may remove or comment these lines in production)
if __name__ == "__main__":
    # Test with a known laurate ester SMILES (e.g., 1-lauroyl-sn-glycerol)
    test_smiles = "CCCCCCCCCCCC(=O)OC[C@@H](O)CO"
    result, reason = is_dodecanoate_ester(test_smiles)
    print(result, reason)
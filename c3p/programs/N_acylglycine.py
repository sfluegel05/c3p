"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine 
Definition: An N-acyl-amino acid in which the amino acid is glycine.
This code checks if a given molecule (in SMILES) contains the substructure:
    R-C(=O)-N-CH2-C(=O)O
which is indicative of N-acylglycine.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    
    An N-acylglycine is an N-acyl amino acid where the amino acid is glycine.
    This function searches for the substructure:
        [CX3](=O)[NX3]CC(=O)[O]
    which represents an acyl group (R-C(=O)) attached to a nitrogen;
    that same nitrogen then bonds to a CH2 group which carries a carboxylic acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains the N-acylglycine substructure, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for N-acylglycine.
    # The pattern [CX3](=O) matches a carbonyl carbon in an acyl group.
    # [NX3] matches a trivalent nitrogen (amide nitrogen).
    # CC(=O)[O] matches a CH2 group attached to a carboxyl (either -COOH or -COO-).
    pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]CC(=O)[O]")
    if pattern is None:
        return False, "Failed to create SMARTS pattern for N-acylglycine"
    
    # Check if the molecule has the N-acylglycine substructure
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains the N-acylglycine substructure (R-C(=O)-N-CH2-C(=O)O)"
    else:
        return False, "N-acylglycine substructure not found"

# Example usage (uncomment for testing):
# test_smiles = "CC(=O)NCC(O)=O"  # N-acetylglycine
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)
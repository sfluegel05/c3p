"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid containing at least one C=C double bond.
A fatty acid is identified by the presence of a carboxylic acid group and a long aliphatic chain.
For this classifier, we require the molecule to have a carboxyl group (COOH or its deprotonated form)
and at least one carbon–carbon double bond.
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    
    A fatty acid is characterized by having a carboxylic acid moiety.
    It is olefinic if it further contains at least one C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an olefinic fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for a carboxylic acid group.
    # We account for both protonated (COOH) and deprotonated (COO-) forms.
    ca_protonated = Chem.MolFromSmarts("C(=O)[O;H1]")  # typical carboxylic acid
    ca_deprotonated = Chem.MolFromSmarts("C(=O)[O-]")   # deprotonated form
    
    # Check for presence of at least one carboxyl group.
    if not (mol.HasSubstructMatch(ca_protonated) or mol.HasSubstructMatch(ca_deprotonated)):
        return False, "No carboxylic acid group found"
    
    # Define SMARTS for a carbon-carbon double bond.
    double_bond = Chem.MolFromSmarts("[C]=[C]")
    
    # Check for at least one C=C double bond.
    if not mol.HasSubstructMatch(double_bond):
        return False, "No carbon–carbon double bond (C=C) found"
    
    return True, "Contains a carboxylic acid group and at least one C=C double bond, classifying it as an olefinic fatty acid"

# Example usage (you can remove or comment these lines when integrating this function into a larger codebase).
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCC\\C=C/CCCCCCCC(O)=O"  # gadoleic acid example
    result, reason = is_olefinic_fatty_acid(test_smiles)
    print(result, reason)
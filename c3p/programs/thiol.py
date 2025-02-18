"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH) is attached 
to a carbon atom (aliphatic or aromatic).
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is defined as an organosulfur compound containing an -SH group attached to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a thiol group, False otherwise
        str: Reason for the classification result
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydrogen on sulfur is properly represented
    mol = Chem.AddHs(mol)

    # Define the SMARTS pattern for a thiol group: a sulfur with one hydrogen attached to a carbon.
    thiol_pattern = Chem.MolFromSmarts("[#6]-[S;H1]")
    
    # Check if the molecule matches the thiol pattern
    if mol.HasSubstructMatch(thiol_pattern):
        return True, "Molecule contains a thiol group (-SH) attached to a carbon atom"
    else:
        return False, "No thiol group (-SH attached to a carbon) found in the molecule"
        
# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example: 2-Methoxybenzenethiol (belongs to the thiol class)
    test_smiles = "SC=1C(OC)=CC=CC1"
    result, reason = is_thiol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)
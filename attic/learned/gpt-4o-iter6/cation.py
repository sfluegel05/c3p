"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    A cation is defined as a monoatomic or polyatomic species having one or more elementary charges of the proton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a positive charge on any atom
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            return True, f"Atom with positive charge found: {atom.GetSymbol()} with charge {atom.GetFormalCharge()}"

    return False, "No positive charge found on any atom"

# Example usage with some SMILES strings for cations
# Note: Actual usage for classification would typically involve loops over datasets
# or interactive/user-provided simulations
print(is_cation("[H]C(=C1C=CN(CCC[N+](C)(C)CCC[N+](C)(C)CCCN2C=CC(C=C2)=C([H])c2oc3ccccc3[n+]2C)C=C1)c1oc2ccccc2[n+]1C"))  # Should output: True
print(is_cation("CCCCC=C"))  # Should output: False
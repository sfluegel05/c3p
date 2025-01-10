"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol has a sulfanyl group (-SH) attached directly to an alkyl group (carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Set up the pattern to find sulfur bound to carbon
    sulfanyl_group = Chem.MolFromSmarts("[S][C]")
    if not mol.HasSubstructMatch(sulfanyl_group):
        return False, "No thiol (-SH) group directly attached to an alkyl carbon found"
    
    # Check for the presence of hydrogen, ensure the sulfur can be a thiol SH 
    atom_with_hydrogen = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16 and atom.GetTotalNumHs() > 0]
    if not atom_with_hydrogen:
        return False, "Sulfur not bonded with a hydrogen atom, not a thiol"

    return True, "Contains alkanethiol pattern with -SH group attached to carbon"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol
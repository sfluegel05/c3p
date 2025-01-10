"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound where a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: Whether the molecule is an alkanethiol
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for thiol group connected to sp3 carbon
    thiol_pattern = Chem.MolFromSmarts("[CX4][SX2H1]")
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No alkanethiol with thiol group (-SH) attached to an alkyl group found"
    
    # Filter out remaining by ensuring sulfur's only single bonds are to hydrogen and carbon
    sulfur_with_hydrogen = [
        atom for atom in mol.GetAtoms() 
        if atom.GetAtomicNum() == 16 and 
           atom.GetTotalNumHs() == 1 and 
           sum(1 for bond in atom.GetBonds() if bond.GetBondTypeAsDouble() == 1) == 2
    ]
    if not sulfur_with_hydrogen:
        return False, "Sulfur is not specifically bonded in sulfanyl (-SH) fashion"

    return True, "Contains alkanethiol pattern with -SH group specifically bound to an alkyl group"

# Example usage
print(is_alkanethiol("SCC(CC)C"))  # Example SMILES for 2-Methyl-1-butanethiol
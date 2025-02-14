"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as having a sulfanyl group (-SH) attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for thiol group attached to an alkyl group
    # The sulfur should be connected to exactly one hydrogen and one carbon
    thiol_pattern = Chem.MolFromSmarts("[SX2H1][#6]")
    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol group (-SH) attached to an alkyl group found"

    # Confirm that the sulfur is not connected to anything beyond the alkyl group
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            connected_atoms = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            # Check if sulfur only bonds with hydrogen and carbon(s)
            if not all(x in ['H', 'C'] for x in connected_atoms):
                return False, f"Sulfur atom bonded with non-alkyl group components: {connected_atoms}"

    return True, "Contains a sulfanyl group, -SH, attached to an alkyl group"
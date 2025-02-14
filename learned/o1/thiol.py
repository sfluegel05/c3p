"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29232 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens (important for accurate hydrogen count)
    mol = Chem.AddHs(mol)

    # Define thiol SMARTS pattern: sulfur with one hydrogen bonded to carbon
    thiol_pattern = Chem.MolFromSmarts("[#16X2H1]-[C]")

    if thiol_pattern is None:
        return False, "Invalid SMARTS pattern for thiol"

    # Search for thiol group
    matches = mol.GetSubstructMatches(thiol_pattern)
    if matches:
        for match in matches:
            s_idx = match[0]  # Index of the sulfur atom
            s_atom = mol.GetAtomWithIdx(s_idx)
            # Verify that sulfur is only bonded to hydrogen and carbon atoms
            neighbor_elements = [neighbor.GetAtomicNum() for neighbor in s_atom.GetNeighbors()]
            if set(neighbor_elements).issubset({1, 6}):
                return True, "Contains a thiol group (-SH) attached to a carbon atom"

        # If no valid thiol is found after checking neighbors
        return False, "No valid thiol group (-SH) attached to a carbon atom found"

    else:
        return False, "No thiol group (-SH) attached to a carbon atom found"
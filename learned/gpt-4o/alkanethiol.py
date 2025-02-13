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

    # Enhanced SMARTS pattern for thiol groups to more accurately target simple alkanethiols
    thiol_pattern = Chem.MolFromSmarts("[#6][SX2H1]")  # -SH must be attached to an sp3 carbon

    if not mol.HasSubstructMatch(thiol_pattern):
        return False, "No thiol group (-SH) directly attached to an alkyl group (sp3 carbon) found"

    # Validate the sulfur atom is not part of larger, non-alkanethiol structures
    for match in mol.GetSubstructMatches(thiol_pattern):
        s_atom = mol.GetAtomWithIdx(match[1])  # Index of sulfur atom in the match
        connected_atoms = [nbr.GetAtomicNum() for nbr in s_atom.GetNeighbors()]

        # Check that sulfur is only attached to carbon (6) and hydrogen (1) atoms
        if any(atom_num not in (1, 6) for atom_num in connected_atoms):
            return False, "Sulfur atom bonded with non-alkyl group components or complex structures"

    return True, "Contains a sulfanyl group, -SH, correctly attached to an alkyl chain"
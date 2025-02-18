"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) indicating if it's an aralkylamine and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amine (N) connected to carbon.
    amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No primary, secondary, or tertiary amine group found"

    # Look for aromatic groups directly bonded to an alkyl carbon and connected to an amine.
    aralkylamine_pattern = Chem.MolFromSmarts("[CX4][a]")  # Alkyl carbon attached directly to an aromatic
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "No alkyl group directly substituted by aromatic group found"

    # Now ensure connection pattern like N-C-C (to ensure we consider the broad range of alkyl connections)
    extended_pattern = Chem.MolFromSmarts("[NX3][CX4][a]")  # N connected to alkyl, then aromatic
    if mol.HasSubstructMatch(extended_pattern):
        return True, "Molecule is an aralkylamine with aromatic substitution on the alkyl side chain"

    return False, "Structure does not match an aralkylamine pattern"
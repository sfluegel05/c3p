"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an aralkylamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for aromatic, alkyl, and amine groups
    aromatic_pattern = Chem.MolFromSmarts("a")  # Aromatic ring
    amine_pattern = Chem.MolFromSmarts("[NX3;!H0]")  # Primary or secondary amine
    aralkyl_linkage_pattern = Chem.MolFromSmarts("a-[CX4]-*!@[NX3]")  # Aromatic connected to alkyl to amine
    
    # Check for the presence of an aromatic ring
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic group found"

    # Check for the presence of an amine group
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"

    # Verify the linkage of the aromatic group via alkyl chain to amine
    if mol.HasSubstructMatch(aralkyl_linkage_pattern):
        return True, "Contains an alkylamine with aromatic substitution"
    else:
        return False, "Aromatic group not properly attached through an alkyl to amine"

    return False, "Not classified as aralkylamine"
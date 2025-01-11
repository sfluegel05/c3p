"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    Omega-hydroxy fatty acids are characterized by a terminal hydroxyl group (ω-position)
    and a carboxylic acid group at the other end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Match omega-terminal hydroxy group
    omega_hydroxy_pattern = Chem.MolFromSmarts("[CH2]CO")
    omega_hydroxy_matches = mol.GetSubstructMatches(omega_hydroxy_pattern)
    if not omega_hydroxy_matches:
        return False, "No omega-terminal hydroxy group found"

    # Match carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Ensure the hydroxy is truly at the omega position in the longest chain
    for match in omega_hydroxy_matches:
        oxy_atom = mol.GetAtomWithIdx(match[1])
        is_omega = all(nb.GetIdx() == match[0] or nb.GetAtomicNum() != 6 for nb in oxy_atom.GetNeighbors())
        if is_omega:
            break
    else:
        return False, "Hydroxy group is not at the omega position"

    # Chain length criteria, allow shorter as examples suggest
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) < 6:
        return False, f"Carbon chain too short, found {len(carbon_atoms)} carbons"

    return True, "Contains omega-terminal hydroxy and carboxylic acid groups with sufficient chain length"
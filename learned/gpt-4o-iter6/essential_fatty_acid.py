"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group (-C(=O)O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Count total number of carbons to ensure minimum chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Insufficient aliphatic chain length (found {carbon_count} carbons), need at least 16"

    # Identify cis double bonds (Z stereochemistry, indicated as \C=C)
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C\C")
    cis_double_bond_count = len(mol.GetSubstructMatches(cis_double_bond_pattern))
    if cis_double_bond_count < 2:
        return False, f"Insufficient cis double bonds (found {cis_double_bond_count}), need at least 2 for polyunsaturation"

    # Check if the molecule has a mostly linear structure, typical of fatty acids
    linearness_threshold_ratio = 0.8  # Expect at least 80% of chain length is aliphatic
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable / carbon_count < linearness_threshold_ratio:
        return False, "Molecule appears not to have a sufficiently linear aliphatic chain"

    return True, "Contains key characteristics of essential fatty acid: carboxylic acid group, multiple cis double bonds, and sufficiently long linear aliphatic chain"
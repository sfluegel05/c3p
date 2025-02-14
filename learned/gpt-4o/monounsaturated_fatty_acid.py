"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain, with singly bonded carbon atoms otherwise.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count double and triple bonds
    doub_triple_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() in [2.0, 3.0])
    
    if doub_triple_bond_count != 1:
        return False, f"Found {doub_triple_bond_count} double/triple bonds, require exactly one"

    # Check that the rest of the chain involves single bonds
    all_single_bonded = all(bond.GetBondTypeAsDouble() == 1.0 or bond.GetBondTypeAsDouble() in [2.0, 3.0] if bond.IsInRing() else False for bond in mol.GetBonds())
    
    if not all_single_bonded:
        return False, "Not all other bonds are single carbon-carbon bonds"

    return True, "Contains one double/triple bond and all other carbon-carbon bonds are single, with a carboxylic acid group"
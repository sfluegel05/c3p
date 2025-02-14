"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    A MUFA has a single carbon chain with one double or triple bond and a terminal carboxylic acid.

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

    # Check for the presence of a terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Identify all carbon atoms part of the main chain (with carboxylic acid)
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Possible sites for attachments to form a chain
    main_chain_carbons = set()

    # Check for main chain including the carboxylic acid group
    for match in carboxylic_acid_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                main_chain_carbons.add(idx)

    # Look for unsaturation (double/triple bonds) within main carbon chain
    unsaturation = 0
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() in [2.0, 3.0]:  # Double or triple bonds
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in main_chain_carbons and end_idx in main_chain_carbons:
                unsaturation += 1

    if unsaturation == 1:
        return True, "Molecule is a monounsaturated fatty acid with one double/triple bond in the main carbon chain"

    return False, f"Found {unsaturation} unsaturations in the potential main chain, require exactly one"
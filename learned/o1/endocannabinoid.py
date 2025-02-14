"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are endogenous ligands that activate cannabinoid receptors,
    typically consisting of a long-chain polyunsaturated fatty acid connected
    via an amide or ester linkage to an ethanolamine or glycerol moiety.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define SMARTS patterns
    # Long-chain fatty acid (at least 15 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")  # 15 carbons

    # At least two double bonds in the chain
    double_bond_pattern = Chem.MolFromSmarts("C=CC=CC=C")  # Multiple conjugated double bonds

    # Amide linkage to ethanolamine (-C(=O)NCCO)
    ethanolamine_amide_pattern = Chem.MolFromSmarts("C(=O)NCCO")

    # Ester linkage to glycerol (-C(=O)OCC(O)CO)
    glycerol_ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]C(O)CO")

    # Ether linkage to glycerol (-COCC(O)CO)
    glycerol_ether_pattern = Chem.MolFromSmarts("COC[C@H](O)CO")

    # Check for long-chain fatty acid moiety
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No long-chain fatty acid moiety found."

    # Check for multiple double bonds in the fatty acid chain
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No polyunsaturated fatty acid chain found."

    # Check for amide linkage to ethanolamine
    if mol.HasSubstructMatch(ethanolamine_amide_pattern):
        return True, "Contains amide linkage to ethanolamine and polyunsaturated fatty acid chain."

    # Check for ester linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ester_pattern):
        return True, "Contains ester linkage to glycerol and polyunsaturated fatty acid chain."

    # Check for ether linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ether_pattern):
        return True, "Contains ether linkage to glycerol and polyunsaturated fatty acid chain."

    return False, "Molecule does not match common endocannabinoid structures."
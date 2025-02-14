"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are endogenous ligands that activate cannabinoid receptors,
    typically consisting of long-chain polyunsaturated fatty acids (usually around 20 carbons)
    connected via an amide linkage to ethanolamine or an ester/ether linkage to glycerol.

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

    # Define SMARTS patterns for identification
    # Ethanolamine amide linkage pattern
    ethanolamine_amide_pattern = Chem.MolFromSmarts("""
        [$([NX3;$([H0,$(H1)])])]-C(=O)-[C;H2]-[C;H2]-[OX2H]
    """)

    # Glycerol ester linkage pattern
    glycerol_ester_pattern = Chem.MolFromSmarts("""
        [$([CX3](=O)[OX2])]-[C;H]([OX2H])-[C;H]([OX2H])-[OX2H]
    """)

    # Glycerol ether linkage pattern
    glycerol_ether_pattern = Chem.MolFromSmarts("""
        [$([OX2])] -[C;H]([OX2H])-[C;H]([OX2H])-[OX2H]
    """)

    # Long-chain polyunsaturated fatty acid pattern (at least 16 carbons, multiple double bonds)
    fatty_acid_chain_pattern = Chem.MolFromSmarts("""
        [CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[*]
    """)

    # Check for fatty acid chain connected to the linkage
    matches_fatty_acid_chain = mol.HasSubstructMatch(fatty_acid_chain_pattern)
    if not matches_fatty_acid_chain:
        return False, "No long-chain fatty acid chain found (at least 16 carbons)."

    # Count the number of double bonds in the fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    double_bond_count = len(double_bond_matches)
    if double_bond_count < 2:
        return False, "Insufficient number of double bonds in fatty acid chain (less than 2)."

    # Check for ethanolamine amide linkage with fatty acid chain
    substructures = mol.GetSubstructMatches(ethanolamine_amide_pattern)
    for match in substructures:
        ethanolamine_amide = Chem.PathToSubmol(mol, match)
        # Check connectivity to fatty acid chain
        if ethanolamine_amide.HasSubstructMatch(fatty_acid_chain_pattern):
            return True, "Contains amide linkage to ethanolamine and polyunsaturated fatty acid chain."

    # Check for glycerol ester linkage with fatty acid chain
    substructures = mol.GetSubstructMatches(glycerol_ester_pattern)
    for match in substructures:
        glycerol_ester = Chem.PathToSubmol(mol, match)
        # Check connectivity to fatty acid chain
        if glycerol_ester.HasSubstructMatch(fatty_acid_chain_pattern):
            return True, "Contains ester linkage to glycerol and polyunsaturated fatty acid chain."

    # Check for glycerol ether linkage with fatty acid chain
    substructures = mol.GetSubstructMatches(glycerol_ether_pattern)
    for match in substructures:
        glycerol_ether = Chem.PathToSubmol(mol, match)
        # Check connectivity to fatty acid chain
        if glycerol_ether.HasSubstructMatch(fatty_acid_chain_pattern):
            return True, "Contains ether linkage to glycerol and polyunsaturated fatty acid chain."

    return False, "Molecule does not match endocannabinoid structures."
"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids are endogenous ligands that activate cannabinoid receptors,
    typically consisting of long-chain polyunsaturated fatty acids connected
    via an amide linkage to ethanolamine or an ester linkage to glycerol.

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

    # Long-chain fatty acid (16-22 carbons) with multiple double bonds
    fatty_acid_pattern = Chem.MolFromSmarts("""
        [#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]-[CH2]-[#6]
        """)  # Chain of at least 16 carbons
    # Pattern for double bonds in the fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("C=CC")  # Double bonds in chain

    # Amide linkage to ethanolamine (-C(=O)NCCO)
    ethanolamine_amide_pattern = Chem.MolFromSmarts("C(=O)NCCO")

    # Ester linkage to glycerol (-C(=O)OCC(O)CO)
    glycerol_ester_pattern = Chem.MolFromSmarts("C(=O)OCC(O)CO")

    # Ether linkage to glycerol (-COCC(O)CO)
    glycerol_ether_pattern = Chem.MolFromSmarts("COCC(O)CO")

    # General glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")

    # Check for long-chain fatty acid moiety
    mol_Hs = Chem.AddHs(mol)
    matches = mol_Hs.GetSubstructMatches(fatty_acid_pattern)
    if not matches:
        return False, "No long-chain fatty acid moiety found."

    # Check for multiple double bonds in the fatty acid chain
    num_double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    if num_double_bonds < 2:
        return False, "No polyunsaturated fatty acid chain found."

    # Check for amide linkage to ethanolamine
    if mol.HasSubstructMatch(ethanolamine_amide_pattern):
        return True, "Contains amide linkage to ethanolamine and polyunsaturated fatty acid chain."

    # Check for ester linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ester_pattern) and mol.HasSubstructMatch(glycerol_pattern):
        return True, "Contains ester linkage to glycerol and polyunsaturated fatty acid chain."

    # Check for ether linkage to glycerol
    if mol.HasSubstructMatch(glycerol_ether_pattern):
        return True, "Contains ether linkage to glycerol and polyunsaturated fatty acid chain."

    return False, "Molecule does not match common endocannabinoid structures."
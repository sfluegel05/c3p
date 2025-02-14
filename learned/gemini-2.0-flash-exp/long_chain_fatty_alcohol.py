"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol (C13-C22) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty alcohol, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a single terminal alcohol group
    # The [CH2X4,CHX4] specifies the carbon must be sp3 hybridized to prevent aromatic or other types of OH groups
    terminal_alcohol_pattern = Chem.MolFromSmarts("[CH2X4,CHX4][OX2H]")
    alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)
    if len(alcohol_matches) != 1:
        return False, f"Molecule has {len(alcohol_matches)} terminal alcohol groups, should have exactly 1."


    # Define SMARTS for a linear alkyl chain attached to the terminal alcohol, with 12 to 21 carbons.
    # This covers C13 to C22.
    # The ~[CH2X4] means that the chain can contain any number of CH2 groups
    # The {11,20} denotes the range for the number of repetition of the preceding atom pattern.
    chain_smarts = "[CH2X4,CHX4]([OX2H])~[CH2X4]{11,20}[CH3X4]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)

    if not mol.HasSubstructMatch(chain_pattern):
      return False, f"Molecule does not have a linear alkyl chain (C13-C22) with a terminal alcohol."

    return True, "Molecule is a long-chain fatty alcohol (C13-C22) with one alcohol group at the end of the chain."
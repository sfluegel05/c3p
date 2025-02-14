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

    # SMARTS pattern for a terminal alcohol attached to a carbon chain
    # The [CX4] indicates a carbon with 4 bonds (sp3 hybridized)
    terminal_alcohol_pattern = Chem.MolFromSmarts("[CH2X4,CHX4][OX2H]")
    
    # Check for the substructure using a loop to ensure single alcohol
    matches = mol.GetSubstructMatches(terminal_alcohol_pattern)
    if len(matches) != 1:
      return False, f"Molecule has {len(matches)} terminal alcohol groups, it should have exactly 1"

    # Count the number of carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check if the chain length is within the range of C13-C22
    if carbon_count < 13 or carbon_count > 22:
        return False, f"Molecule has {carbon_count} carbons, should be between 13 and 22"

    # Check for long alkyl chain
    chain_pattern = Chem.MolFromSmarts("C[C,CH2,CH]") # This pattern defines a chain of C, CH and CH2.
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 10:
      return False, "The chain is too short to qualify as a fatty alcohol."

    return True, "Molecule is a long-chain fatty alcohol (C13-C22) with one alcohol group at the end of the chain."
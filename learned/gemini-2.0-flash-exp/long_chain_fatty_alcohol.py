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

    # Check for at least one alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[CX4,CX3][OX2H]")
    if alcohol_pattern is None:
        return False, "Could not parse alcohol SMARTS pattern"

    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) == 0:
        return False, "Molecule does not have an alcohol group"

    # Identify a long chain (at least 12 carbons) using SMARTS
    # Use a more flexible pattern to identify a chain within the molecule
    chain_smarts = "[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    if chain_pattern is None:
        return False, "Could not parse chain SMARTS pattern"

    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Molecule does not contain a chain with at least 12 carbons."


    # Iterate through all atoms and check for carbon chain length and ensure that it is
    # between 13 and 22 and that it contains a hydroxyl group.
    
    chain_carbon_count = 0
    has_oh = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # Find all non-ring atoms
            if not atom.IsInRing():
                chain_carbon_count += 1
        if atom.GetAtomicNum() == 8 and atom.HasProp('molAtomMapNumber') == False:
            has_oh = True

    if not has_oh:
            return False, "Molecule has no hydroxyl group"


    if chain_carbon_count < 13:
          return False, f"Molecule has a carbon chain less than 13 carbons ({chain_carbon_count})"
    if chain_carbon_count > 22:
         return False, f"Molecule has a carbon chain more than 22 carbons ({chain_carbon_count})"

    return True, f"Molecule is a long-chain fatty alcohol (C13-C22) with {chain_carbon_count} carbon atoms."
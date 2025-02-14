"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related
    to cholestan-3-ol, with possible additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broadened sterol core pattern (includes cholestane framework and flexibility for side chains)
    sterol_core_smarts = "[#6]1([#6][#6][#6]2)[#6][#6][#6]3[#6](C)([#6][#6][#6]4)[#6]([OH])[#6]4[#6][#6][#6]3[#6]2C1"
    sterol_core_pattern = Chem.MolFromSmarts(sterol_core_smarts)
    if not mol.HasSubstructMatch(sterol_core_pattern):
        return False, "Does not contain a sterol core structure"

    # Validate presence of a 3-hydroxy group
    # The hydroxy group should appear at a specific position related to the sterol core
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[#6]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have a 3-hydroxy group"

    # Check alignment with a broader sterol skeleton potentially including additional chains
    additional_chain_smarts = "[C@H][C]"
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts(additional_chain_smarts))) < 3:
        return True, "Recognized as a sterol with 3-hydroxy group and adaptable skeleton"

    return False, "Failed to match sterol's unique skeletal and functional group requirements"
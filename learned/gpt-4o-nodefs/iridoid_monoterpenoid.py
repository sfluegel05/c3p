"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically contains a cyclopentanopyran framework and may
    include characteristic functional groups such as lactones or carboxylic esters.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the iridoid monoterpenoid class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS patterns for cyclopentanopyran frameworks and related structures
    cyclopentanopyran_smarts_variants = [
        "C1C(O)C=CC2C1CC=C2",  # Cyclopentanopyran with different attachments
        "C1CC2C(C1)OC=C2",  # Traditional cyclopentanopyran core
        "C1C=CC(C1)C=O",  # Open chain forms leading to cyclic structures
        "[O]1C=CC2C(O)C(C)CC2C1(O)  # More involved shard structures with oxygen bridges
    ]

    # Check for common functional groups associated with iridoid monoterpenoids
    functional_groups_smarts = [
        "[CX3](=O)[OX2H1]",  # Carboxylic acid
        "[CX3](=O)[O][CX3](=O)[OX2H1]",  # Ester linkage
        "C[C@H]1CC[C@@H](O)[C@H](O)C(=O)O1"  # Lactone ring
    ]

    # Attempt to match any of the SMARTS patterns to the molecule
    for cp_smarts in cyclopentanopyran_smarts_variants:
        pattern = Chem.MolFromSmarts(cp_smarts)
        if mol.HasSubstructMatch(pattern):
            # Further validate the presence of functional groups
            for fg_smarts in functional_groups_smarts:
                fg_pattern = Chem.MolFromSmarts(fg_smarts)
                if mol.HasSubstructMatch(fg_pattern):
                    return True, "Compound contains a cyclopentanopyran-related framework with characteristic functional groups"

    return False, "No cyclopentanopyran or related framework with characteristic functional groups found; does not match typical iridoid monoterpenoid structure"
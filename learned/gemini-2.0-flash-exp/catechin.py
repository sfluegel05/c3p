"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ols, having the core structure of two benzene rings
    linked by a pyran ring and a hydroxyl group at the 3-position of the pyran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catechin, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core flavan-3-ol substructure (generic, no specific substitutions).
    flavan_3ol_core = Chem.MolFromSmarts("[CH2][CH](O)[CH2]1[c]2[c]([c][c][c][c]2)[O]1~[c]3[c][c][c][c][c]3")
    if not mol.HasSubstructMatch(flavan_3ol_core):
        return False, "Core flavan-3-ol structure not found."

    # Exclude flavanones that have a ketone group at position 4.
    flavanone_core = Chem.MolFromSmarts("[CH2][C](=[O])[CH2]1[c]2[c]([c][c][c][c]2)[O]1")
    if mol.HasSubstructMatch(flavanone_core):
      return False, "Molecule is a flavanone, not a catechin"

    return True, "Molecule contains the core flavan-3-ol structure."
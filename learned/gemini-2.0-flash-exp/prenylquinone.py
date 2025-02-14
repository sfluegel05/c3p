"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a quinone core (a conjugated ring system with two carbonyls in a 1,4 position)
    quinone_pattern = Chem.MolFromSmarts("[C;$(C=O)]1[C,c](=[C,c][C,c](=[C,c]1)[C;$(C=O)])")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone core found"

    # Check for at least two isoprenoid units (C-C=C-C) attached to a carbon
    #Prenyl group is defined as a carbon chain containing a double bond (C=C)
    prenyl_pattern = Chem.MolFromSmarts("[CX4]!@[CH2]C=C[CX4]")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if len(prenyl_matches) < 2:
        return False, "Less than two polyprenyl side chains found"


    return True, "Contains a quinone core and a polyprenyl side chain"
"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is any aliphatic sulfide in which at least one of the organyl groups
    attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide SMARTS pattern
    # [CH3]: methyl group
    # [S;D2]: sulfur atom with degree 2 (connected to two atoms)
    # [#6]: any carbon atom
    methyl_sulfide_pattern = Chem.MolFromSmarts("[CH3]-[S;D2]-[#6]")
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains methyl sulfide group"
    else:
        return False, "Does not contain methyl sulfide group"
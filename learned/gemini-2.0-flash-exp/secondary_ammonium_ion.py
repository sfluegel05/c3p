"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is a nitrogen atom with a positive charge and two directly attached
    carbon containing groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for positively charged nitrogen with 2 directly attached carbon groups
    # Using SMARTS pattern: [NX3+] ([NX3+] = positively charged nitrogen, X3 = 3 bonds).
    # And we need to specify that these bonds are only to carbon atoms using [CX4,CX3], 
    # where C = carbon, X4 = 4 bonds (sp3 hybridized) or X3 = 3 bonds (sp2 hybridized).
    ammonium_pattern = Chem.MolFromSmarts("[NX3+]([CX4,CX3])[CX4,CX3]")

    if mol.HasSubstructMatch(ammonium_pattern):
            return True, "Contains a positively charged nitrogen with 2 directly attached carbon groups."
    else:
        return False, "Does not contain a positively charged nitrogen with 2 directly attached carbon groups."
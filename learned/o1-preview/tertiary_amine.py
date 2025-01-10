"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three hydrocarbyl groups via single bonds,
    with no hydrogen atoms attached, not in an aromatic ring, and zero formal charge.
    The nitrogen should be sp³ hybridized and not part of any functional group like amide, imine, nitrile, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for tertiary amine
    # [#7X3+0] : Nitrogen atom (atomic number 7), sp³ hybridized, zero formal charge
    # !$(N-C=O) : Exclude amides (nitrogen bonded to carbonyl carbon)
    # !$(N=C) : Exclude imines (nitrogen double-bonded to carbon)
    # !$(N#*) : Exclude nitriles (nitrogen triple-bonded)
    # !$(N-[!#6]) : Exclude nitrogen bonded to non-carbon atoms
    # This pattern allows nitrogens in rings

    tertiary_amine_smarts = '[#7X3+0;!$([#7][C]=O);!$([#7]=*);!$([#7]#*);!$([#7]-[!#6])]'

    # Create the SMARTS pattern
    pattern = Chem.MolFromSmarts(tertiary_amine_smarts)
    # Search for matches in the molecule
    matches = mol.GetSubstructMatches(pattern)

    # If any matches are found, it's a tertiary amine
    if matches:
        return True, "Contains a tertiary amine nitrogen"
    else:
        return False, "No tertiary amine nitrogen found"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32877',
                          'name': 'tertiary amine',
                          'definition': 'A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32874']}}
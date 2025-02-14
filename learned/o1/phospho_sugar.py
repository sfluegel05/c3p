"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for phosphate ester group attached to an alcohol
    phosphate_ester_pattern = Chem.MolFromSmarts("[$([OX2H0][P](=O)([OX1])[OX1])]")

    # Search for phosphate ester groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester_pattern)
    if not phosphate_matches:
        return False, "No phosphate ester group attached to an alcohol found"

    # Define a SMARTS pattern for a monosaccharide (3 to 7 carbon backbone with hydroxyl groups)
    sugar_pattern = Chem.MolFromSmarts("""
        [$([C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)),
         $([C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)),
         $([C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)),
         $([C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)),
         $([C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O)[C;H1,H2](O))]
    """)

    # Search for sugar moiety
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No monosaccharide moiety found"

    # Check if the phosphate ester is connected to the sugar
    phosphate_atoms = set()
    for match in phosphate_matches:
        phosphate_atoms.update(match)

    sugar_atoms = set()
    for match in sugar_matches:
        sugar_atoms.update(match)

    # Find common atoms between phosphate ester group and sugar moiety
    common_atoms = phosphate_atoms.intersection(sugar_atoms)
    if not common_atoms:
        return False, "Phosphate ester group is not attached to the sugar moiety"

    # If passes all checks
    return True, "Contains monosaccharide with an alcoholic hydroxy group esterified with phosphoric acid"
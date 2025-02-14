"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for sphingoid backbone
    # This pattern captures:
    # - A long aliphatic chain (at least 12 carbons)
    # - An amino group at C2 (can be primary amine or amide)
    # - Hydroxyl groups at C1 and C3 positions
    # - Allows for variations like double bonds and additional hydroxyl groups
    sphingoid_pattern = Chem.MolFromSmarts("""
    [#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]
    """)
    # Simplify the pattern to focus on key features
    sphingoid_pattern = Chem.MolFromSmarts("""
    [#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](-[O])[#6]-[#6]-[#6]-[#7X3]-[#6](-[O])-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]
    """)

    if sphingoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Find matches for the sphingoid backbone
    matches = mol.GetSubstructMatches(sphingoid_pattern)

    if not matches:
        return False, "Molecule does not match sphingoid backbone"

    # Check each match for chain length and functional groups
    for match in matches:
        atom_indices = list(match)
        submol = Chem.PathToSubmol(mol, atom_indices)

        # Check for chain length (at least 12 carbons)
        carbon_atoms = [atom for atom in submol.GetAtoms() if atom.GetAtomicNum() == 6]
        if len(carbon_atoms) < 12:
            continue  # Chain too short to be sphingoid

        # Check for amino group at C2 (primary amine or amide)
        nitrogen_atoms = [atom for atom in submol.GetAtoms() if atom.GetAtomicNum() == 7]
        if not nitrogen_atoms:
            continue  # No nitrogen found

        # Check for hydroxyl groups at C1 and C3 positions
        oxygen_atoms = [atom for atom in submol.GetAtoms() if atom.GetAtomicNum() == 8]
        if len(oxygen_atoms) < 2:
            continue  # Not enough hydroxyl groups

        # Additional checks can be added here if necessary

        # If all checks passed
        return True, "Molecule matches sphingoid structural features"

    return False, "Molecule does not match sphingoid structural features"
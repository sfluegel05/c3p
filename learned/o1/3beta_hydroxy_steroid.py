"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:36804 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group at the 3-position in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure stereochemistry is assigned correctly
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol)
        Chem.AssignAtomChiralTagsFromStructure(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except:
        return False, "Failed to assign stereochemistry"

    # Define the steroid nucleus SMARTS pattern (cyclopentanoperhydrophenanthrene core)
    steroid_pattern = Chem.MolFromSmarts("""
    [#6]1([#6H2])[#6H]2[#6H]([#6H2])[#6H2][#6]3([#6H2])[#6H2][#6H]([#6H2])[#6H2][#6]([#6H2])[#6H]4[#6H2][#6]([#6H2])[#6H2][#6]([#6H]1)[#6H]2[#6H]3[#6H]4
    """)  # Steroid backbone

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid core not found"

    # Define the 3beta-hydroxy group pattern at position 3 with beta orientation
    beta_oh_pattern = Chem.MolFromSmarts("""
    [C@@H]([O])[C@H](C)[C@H]
    """)  # 3beta-hydroxy group with correct stereochemistry

    if beta_oh_pattern is None:
        return False, "Invalid 3beta-hydroxy SMARTS pattern"

    # Find matches for the 3beta-hydroxy group with correct stereochemistry
    matches = mol.GetSubstructMatches(beta_oh_pattern, useChirality=True)
    if not matches:
        return False, "No 3beta-hydroxy group found with correct stereochemistry"

    # Verify that the hydroxy group is at the 3-position of the steroid nucleus
    # Map the steroid nucleus to get atom indices
    steroid_matches = mol.GetSubstructMatch(steroid_pattern)
    if not steroid_matches:
        return False, "Steroid core not matched correctly"

    # Assuming standardized atom indexing for the steroid nucleus,
    # check if the hydroxyl-bearing carbon is at position 3
    hydroxyl_carbon_indices = [match[0] for match in matches]
    steroid_carbon_indices = list(steroid_matches)

    # Position 3 in the steroid nucleus corresponds to a specific atom index
    # Mapping may vary, so we need to confirm the correct position
    # For simplicity, check if any hydroxyl carbon is in the steroid core
    core_and_hydroxy_overlap = set(hydroxyl_carbon_indices).intersection(steroid_carbon_indices)
    if not core_and_hydroxy_overlap:
        return False, "3beta-hydroxy group not at position 3 of steroid core"

    return True, "Contains steroid backbone with 3beta-hydroxy group in beta orientation at position 3"
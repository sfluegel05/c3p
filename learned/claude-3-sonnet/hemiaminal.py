"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:48243 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic amino compound that has an amino group and a
    hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hemiaminal substructure
    hemiaminal_pattern = Chem.MolFromSmarts("[NX3,NX4+0][CX4][OX2H]")

    # Find matches for hemiaminal pattern
    matches = mol.GetSubstructMatches(hemiaminal_pattern)

    if not matches:
        return False, "No hemiaminal substructure found"

    # Check stereochemistry and connectivity
    for match in matches:
        nitrogen_idx = match[0]
        carbon_idx = match[1]
        oxygen_idx = match[2]

        # Check if the nitrogen and oxygen are connected to the same carbon
        if not mol.GetBondBetweenAtoms(carbon_idx, nitrogen_idx):
            continue
        if not mol.GetBondBetweenAtoms(carbon_idx, oxygen_idx):
            continue

        # Check stereochemistry (if defined)
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        if carbon_atom.HasProp("_CIPCode"):
            carbon_chirality = int(carbon_atom.GetProp("_CIPCode"))
            if carbon_chirality != Chem.rdchem.CHI_TETRAHEDRAL_CCW:
                continue  # Incorrect stereochemistry

        # Check for other functional groups or substituents
        # that may disqualify the structure from being a hemiaminal
        # ...

        # If all checks pass, it's a hemiaminal
        return True, "Contains an amino group and a hydroxy group attached to the same carbon atom"

    return False, "No valid hemiaminal substructure found"
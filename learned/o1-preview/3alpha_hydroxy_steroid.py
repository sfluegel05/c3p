"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a steroid with a hydroxyl group at position 3 in the alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get Murcko scaffold of the molecule
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not get scaffold of molecule"

    # Define steroid scaffold SMARTS pattern (tetracyclic core with fused rings)
    steroid_smarts = '[R][R][R][R]'  # Molecule with at least 4 rings
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)

    # Check if scaffold matches steroid pattern
    if not scaffold.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"

    # Identify hydroxyl groups attached to ring carbons
    hydroxyl_smarts = '[C@H](O)[C]'
    hydroxyl_pattern = Chem.MolFromSmarts(hydroxyl_smarts)
    matches = mol.GetSubstructMatches(hydroxyl_pattern, useChirality=True)

    if not matches:
        return False, "No chiral hydroxyl group found"

    # Check if any hydroxyl group is in the alpha orientation
    for match in matches:
        carbon_idx = match[0]
        oxygen_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        chiral_tag = carbon_atom.GetChiralTag()

        # Check if attached carbon is part of a ring (steroid ring system)
        if not carbon_atom.IsInRing():
            continue

        # In RDKit, chiral tag 'CHI_TETRAHEDRAL_CCW' corresponds to '@' (alpha)
        # 'CHI_TETRAHEDRAL_CW' corresponds to '@@' (beta)
        if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return True, "3alpha-hydroxy group found with correct stereochemistry"
        else:
            return False, "Hydroxyl group has incorrect stereochemistry (not alpha)"

    return False, "No 3alpha-hydroxy group found"
"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35485 quaternary ammonium ion
A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen
have been replaced with univalent (usually organyl) groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for quaternary nitrogen with four univalent substituents
    quat_n_pattern = Chem.MolFromSmarts("[N+]([H,C,c,O,F])([H,C,c,O,F])([H,C,c,O,F])([H,C,c,O,F])")
    quat_n_matches = mol.GetSubstructMatches(quat_n_pattern)

    if quat_n_matches:
        # Check if the molecule contains a phosphocholine group
        phosphocholine_pattern = Chem.MolFromSmarts("OP(=O)(OCC[N+]([H,C,c,O,F])([H,C,c,O,F])([H,C,c,O,F])([H,C,c,O,F]))[O-]")
        is_phosphocholine = mol.HasSubstructMatch(phosphocholine_pattern)

        # Additional checks
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        n_atoms = mol.GetNumAtoms()
        n_bonds = mol.GetNumBonds()

        if is_phosphocholine or (200 < mol_wt < 1000 and 10 < n_atoms < 100 and 10 < n_bonds < 150):
            return True, "Contains quaternary nitrogen atom with four univalent substituents"

    return False, "No quaternary nitrogen atom with four univalent substituents found"
"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are characterized by repeating isoprene units with a terminal alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic isoprene unit pattern (C=C-C-C core)
    isoprene_core_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX4]-[CX4]")
    core_matches = mol.GetSubstructMatches(isoprene_core_pattern)

    if len(core_matches) < 2:
        return False, f"Too few isoprene core units. found {len(core_matches)}"

    # Check for methyl groups (one per isoprene unit). we count CH3 connected to any carbon in the structure
    methyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)

    if len(methyl_matches) < 2:
       return False, "Too few methyl groups found"

    # Check for terminal alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if len(alcohol_matches) < 1:
        return False, "No terminal alcohol group found"

    # Check for chain length (number of carbon atoms) and rotatable bonds
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    if num_carbons < 10:
      return False, "Too few carbons in the chain"
    if n_rotatable < 2:
        return False, "Too few rotatable bonds for a polyprenol"

    # check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for polyprenol"

    return True, "Molecule has isoprene units, a terminal alcohol and sufficient chain length, compatible with polyprenol definition."
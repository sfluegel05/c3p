"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:34974 fatty amide
A monocarboxylic acid amide derived from a fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (-C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_match = mol.GetSubstructMatch(amide_pattern)
    if not amide_match:
        return False, "No amide group found"

    # Look for carbon chain attached to amide
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]")
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern, amide_match[0])
    if not carbon_chain_matches:
        return False, "No carbon chain found attached to amide group"

    # Check if amide is part of a larger ring system
    ring_info = mol.GetRingInfo()
    amide_atom = mol.GetAtomWithIdx(amide_match[0])
    if any(amide_atom.IsInRingSize(size) for size in ring_info.AtomRings()):
        return False, "Amide group is part of a ring system"

    # Check for common functional groups found in lipids/carbohydrates
    lipid_pattern = Chem.MolFromSmarts("[OX2H,OX1H0-]")
    carbohydrate_pattern = Chem.MolFromSmarts("[OX2H][CX4][OX2H]")
    if mol.HasSubstructMatch(lipid_pattern) or mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "Molecule contains functional groups typical of lipids or carbohydrates"

    # Count rotatable bonds to verify chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Carbon chain too short to be a fatty amide"

    # Check molecular weight - fatty amides typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty amide"

    return True, "Contains amide group with a carbon chain (fatty acid)"
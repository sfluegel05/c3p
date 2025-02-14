"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:46661 lactol

A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic
or ketonic carbonyl group. They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for intramolecular hemiacetal pattern
    hemiacetal_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX3](=[OX1])[CX4][CX4][CX4][CX4]1")
    if not mol.HasSubstructMatch(hemiacetal_pattern):
        return False, "No intramolecular hemiacetal group found"

    # Look for lactone pattern (to exclude lactones)
    lactone_pattern = Chem.MolFromSmarts("[OX2]=[CX3]1[OX2][CX4][CX4][CX4][CX4]1")
    if mol.HasSubstructMatch(lactone_pattern):
        return False, "Lactone structure detected, not a lactol"

    # Check for common structural features of lactols
    unsaturated = mol.GetBondBetweenAtoms(*list(hemiacetal_pattern.GetBondWithStereo(-1).GetBeginAtomIdx())[:-1]).GetBondType() == Chem.BondType.DOUBLE
    has_hydroxy = any(atom.GetAtomicNum() == 8 and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1 for atom in mol.GetAtoms())
    ring_size = len(list(hemiacetal_pattern.GetRingInfo().AtomRings()[0]))

    if unsaturated or not has_hydroxy or ring_size < 4 or ring_size > 10:
        return False, "Structure does not match typical lactol features"

    # Check for reasonable molecular properties
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    if mol_weight < 100 or mol_weight > 800 or n_rotatable > 15:
        return False, "Molecular properties outside expected range for lactols"

    return True, "Contains intramolecular hemiacetal group and matches expected structural features of lactols"
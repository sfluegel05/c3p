from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_bile_acid_glycine_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid glycine conjugate (amide of a bile acid with glycine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid glycine conjugate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid core pattern (more flexible pattern)
    steroid_core = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~[#6]~2~[#6]~1')
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Check for glycine amide group (more flexible pattern)
    glycine_amide = Chem.MolFromSmarts('[#6][CX3](=O)[NX3][CH2][CX3](=O)[OX2H,OX1-]')
    if not mol.HasSubstructMatch(glycine_amide):
        return False, "No glycine amide group found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Insufficient number of rings for bile acid structure"

    # Check for typical number of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient number of oxygen atoms for bile acid glycine conjugate"

    # Verify carbon count
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 24:
        return False, f"Carbon count ({carbon_count}) too low for bile acid glycine conjugates"

    # Check molecular weight
    mol_weight = Descriptors.ExactMolWt(mol)
    if not (350 <= mol_weight <= 800):
        return False, f"Molecular weight {mol_weight:.1f} outside typical range for bile acid glycine conjugates"

    # Additional check for bile acid core with glycine
    bile_acid_glycine = Chem.MolFromSmarts('[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~[#6]~2~[#6]~1~[#6]~[#6]~[#6]~[CX3](=O)[NX3][CH2][CX3](=O)[OX2H,OX1-]')
    if not mol.HasSubstructMatch(bile_acid_glycine):
        return False, "Structure does not match bile acid glycine conjugate pattern"

    return True, "Molecule contains bile acid core conjugated with glycine via amide linkage"
# Pr=None
# Recall=0.0
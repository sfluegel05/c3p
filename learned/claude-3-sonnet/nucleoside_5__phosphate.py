"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:33696 nucleoside 5'-phosphate

A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which C-5
of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ribose/deoxyribose backbone
    ribose_patterns = ["[OX2]C1[CH]([OX2])[CH]([OX2])[CH]([OX2])[CH]1",
                       "[OX2]C1[CH]([OX2])[CH]([OX2])[CH]([OX2])C1"]
    has_ribose = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in ribose_patterns)
    if not has_ribose:
        return False, "No ribose/deoxyribose backbone found"

    # Check for pyrimidine/purine base
    base_patterns = ["c1ncnc[nH]1", "c1[nH]cnc1", "n1cnc2c1ncnc2", "c1ncnc2[nH]cnc12"]
    has_base = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in base_patterns)
    if not has_base:
        return False, "No pyrimidine/purine base found"

    # Check for phosphate group(s) attached to C-5 of ribose/deoxyribose
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O[C@H]1[C@@H]([CH]([CH]([CH]1[OX2])O)O)O")
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group attached to C-5 of ribose/deoxyribose"

    # Additional checks
    num_phosphates = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if num_phosphates < 1 or num_phosphates > 4:
        return False, "Number of phosphate groups must be between 1 and 4"

    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300 or mol_weight > 1000:
        return False, "Molecular weight outside typical range for nucleoside phosphates"

    return True, "Molecule satisfies structural requirements for nucleoside 5'-phosphate"
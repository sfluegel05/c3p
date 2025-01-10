"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    A 3-oxo-fatty acyl-CoA(4-) is an acyl-CoA(4-) with a 3-oxo-fatty acyl chain attached via a thioester linkage,
    where the phosphate and diphosphate groups are deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester linkage with 3-oxo-fatty acyl chain
    # Pattern: O=C-C(=O)-S-C (beta-keto group connected via thioester linkage)
    beta_keto_thioester_smarts = Chem.MolFromSmarts('O=C-C(=O)-S-C')
    if not mol.HasSubstructMatch(beta_keto_thioester_smarts):
        return False, "3-oxo-fatty acyl chain with thioester linkage not found"

    # Check for pantetheine unit (simplified)
    # Pattern: NC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP
    pantetheine_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP')
    if not mol.HasSubstructMatch(pantetheine_smarts):
        return False, "Pantetheine unit not found"

    # Check for adenosine diphosphate moiety (simplified)
    # Pattern: n1cnc2c1ncnc2N
    adenosine_smarts = Chem.MolFromSmarts('n1cnc2c1ncnc2N')
    if not mol.HasSubstructMatch(adenosine_smarts):
        return False, "Adenosine moiety not found"

    # Check for diphosphate groups (deprotonated)
    # Pattern: [O-]P(=O)([O-])OP(=O)([O-])[O-]
    diphosphate_smarts = Chem.MolFromSmarts('[O-]P(=O)([O-])OP(=O)([O-])[O-]')
    if not mol.HasSubstructMatch(diphosphate_smarts):
        return False, "Deprotonated diphosphate groups not found"

    # Check for net charge of -4
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Incorrect net charge ({total_charge}), expected -4 for deprotonated phosphate groups"

    # The molecule passes all checks
    return True, "Molecule is a 3-oxo-fatty acyl-CoA(4-) with deprotonated phosphate groups"

__metadata__ = {
    'chemical_class': {
        'name': '3-oxo-fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups of any 3-oxo-fatty acyl-CoA.',
        'examples': [
            '3-oxoisooctadecanoyl-CoA(4-)',
            '3-oxohexacosanoyl-CoA(4-)',
            '3-oxooctanoyl-CoA(4-)',
            # ... (additional examples)
        ]
    }
}
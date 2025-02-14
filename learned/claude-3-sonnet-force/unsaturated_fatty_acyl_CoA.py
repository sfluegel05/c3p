"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: CHEBI:35497 unsaturated fatty acyl-CoA
A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any unsaturated fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA backbone
    coa_pattern = Chem.MolFromSmarts("[N;X3]1([C;X4][N;X3]([C;X3]([C;X3]([C;X3]([C;X3]1)O[P;X4]([O-])([O-])OP([O-])([O-])O[CH2]O[CH2]O[P;X4]([O-])([O-])OP([O-])([O-])OC[C;X4][N;X3]c1[nH]c[nH]c1)[N;X3]c1[nH]c[nH]c1)[N;X3]c1[nH]c[nH]c1)[O-])[C@H]([C@H]([C@@H]([C@H](O)O)O)O)OP([O-])([O-])[O-]")
    if mol.HasSubstructMatch(coa_pattern):
        reason = "Contains CoA backbone"
    else:
        return False, "Missing CoA backbone"

    # Check for unsaturated alkyl chain
    alkyl_pattern = Chem.MolFromSmarts("[CH2]=[CH][CH2][CH2][CH2][CH2][CH2]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    if alkyl_matches:
        reason = f"{reason}, found unsaturated alkyl chain"
    else:
        return False, "No unsaturated alkyl chain found"

    # Check for carboxylate group (fatty acid)
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if carboxylate_matches:
        reason = f"{reason}, and carboxylate group present (fatty acid)"
    else:
        return False, "No carboxylate group found (not a fatty acid)"

    # Check for unsaturated bond stereochemistry
    unsaturated_bond_stereo = any(bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if unsaturated_bond_stereo:
        reason = f"{reason}, with defined stereochemistry"
    else:
        reason = f"{reason}, without defined stereochemistry"

    # Optionally check molecular weight
    # mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # if mol_wt < 500:
    #     return False, "Molecular weight too low for unsaturated fatty acyl-CoA"

    return True, reason
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

    # Check for Coenzyme A moiety (pantetheine + adenosine diphosphate)
    # Using a SMARTS pattern that matches the CoA structure without the terminal thiol group
    coa_smarts = Chem.MolFromSmarts("""
        [NH2][C](=O)[C][C][NH][C](=O)[C@@H]([O])[C]([C])([C])[C][OP](=O)([O-])[O][P](=O)([O-])[O][C][C@@H]1[O][C@@H]([n]2[c][n][c]3[c][n][c]([NH2])[n][c]23)[C@@H]([O])[C@@H]1[OP](=O)([O-])[O-]
    """)
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Find the thioester linkage: C(=O)-S-
    thioester_smarts = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage found
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[1]

        # Get the alpha carbon (next carbon in acyl chain)
        alpha_c = None
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() != sulfur_idx and neighbor.GetAtomicNum() == 6:
                alpha_c = neighbor
                break
        if alpha_c is None:
            continue  # Can't find alpha carbon, move to next match

        # Get the beta carbon (second carbon in acyl chain)
        beta_c = None
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_c_idx and neighbor.GetAtomicNum() == 6:
                beta_c = neighbor
                break
        if beta_c is None:
            continue  # Can't find beta carbon, move to next match

        # Check if the beta carbon has a carbonyl group (C=O)
        has_beta_keto = False
        for neighbor in beta_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and beta_c.GetBondBetween(beta_c, neighbor).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_beta_keto = True
                break
        if not has_beta_keto:
            continue  # Beta carbon does not have keto group, move to next match

        # If we reach here, the molecule has a 3-oxo-fatty acyl chain attached via thioester linkage
        # Now, check for deprotonated phosphate groups (net charge -4)
        total_charge = Chem.GetFormalCharge(mol)
        if total_charge != -4:
            return False, f"Incorrect net charge ({total_charge}), expected -4 for deprotonated phosphate groups"

        # All checks passed
        return True, "Molecule is a 3-oxo-fatty acyl-CoA(4-) with deprotonated phosphate groups"

    # If no thioester linkage with beta-keto group found
    return False, "No 3-oxo-fatty acyl chain attached via thioester linkage found"

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
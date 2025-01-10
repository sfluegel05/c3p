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

    # Coenzyme A moiety SMILES (deprotonated form)
    coa_smiles = 'NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(n2cnc3c(N)ncnc23)C(O)C1OP(=O)([O-])[O-]'
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Invalid CoA SMILES pattern"

    # Check for Coenzyme A moiety
    if not mol.HasSubstructMatch(coa_mol):
        return False, "Coenzyme A moiety not found"

    # Find thioester linkage: C(=O)-S-C (between acyl chain and CoA)
    thioester_smarts = 'C(=O)SC'
    thioester_mol = Chem.MolFromSmarts(thioester_smarts)
    if thioester_mol is None:
        return False, "Invalid thioester SMARTS pattern"

    thioester_matches = mol.GetSubstructMatches(thioester_mol)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage found
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # C in C(=O)
        sulfur_idx = match[2]      # S in C-S-C
        beta_c = None
        gamma_c = None

        # Get the alpha carbon (adjacent to carbonyl carbon in acyl chain)
        alpha_c = None
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        for neighbor in carbonyl_c.GetNeighbors():
            if neighbor.GetIdx() != match[1] and neighbor.GetAtomicNum() == 6:
                alpha_c = neighbor
                break
        if alpha_c is None:
            continue  # Can't find alpha carbon, move to next match

        # Get the beta carbon (next carbon in acyl chain)
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetIdx() != carbonyl_c.GetIdx() and neighbor.GetAtomicNum() == 6:
                beta_c = neighbor
                break
        if beta_c is None:
            continue  # Can't find beta carbon, move to next match

        # Check if the beta carbon has a carbonyl group (C=O) attached (3-oxo group)
        has_beta_keto = False
        for neighbor in beta_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                bond = beta_c.GetBondBetween(beta_c, neighbor)
                if bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
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
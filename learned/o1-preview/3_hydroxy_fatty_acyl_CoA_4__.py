"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    A 3-hydroxy fatty acyl-CoA(4-) is an acyl-CoA molecule where the fatty acyl chain has a hydroxyl group at the 3-position,
    and the molecule may be deprotonated at the phosphate and diphosphate groups, resulting in a net charge of -4 or 0.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Standardize the molecule (e.g., kekulize, add hydrogens)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass

    # Check net charge (allowing -4 or 0)
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4 and total_charge != 0:
        return False, f"Total charge is {total_charge}, expected -4 or 0"

    # Check for Coenzyme A moiety
    coenzyme_a_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)[C](C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](n2cnn3c(N)ncnc23)[C@H](O)[C@H]1OP(O)(O)=O')
    if not mol.HasSubstructMatch(coenzyme_a_smarts):
        return False, "Coenzyme A moiety not found"

    # Check for thioester linkage
    thioester_smarts = Chem.MolFromSmarts('C(=O)SCCN')
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "Thioester linkage not found"

    # Check for 3-hydroxy fatty acyl chain attached to thioester
    # The chain attached to the thioester should have a hydroxyl at the 3-position from the carbonyl carbon
    # Define the pattern for the 3-hydroxy fatty acyl chain
    hydroxy_acyl_chain_smarts = Chem.MolFromSmarts('C(=O)SC[C;X4][C@H](O)[C;X4]')
    if not mol.HasSubstructMatch(hydroxy_acyl_chain_smarts):
        return False, "3-hydroxy fatty acyl chain not found"

    # Ensure that the fatty acyl chain is of sufficient length (e.g., at least 2 carbons beyond the hydroxyl)
    # Count the number of carbons in the fatty acyl chain
    acyl_chain_smarts = Chem.MolFromSmarts('C(=O)SC([#6]*)[#6]*')
    acyl_matches = mol.GetSubstructMatches(acyl_chain_smarts)
    if acyl_matches:
        # Get the fatty acyl chain fragment
        acyl_atom_indices = acyl_matches[0]
        acyl_atom = mol.GetAtomWithIdx(acyl_atom_indices[2])  # The carbon attached to the sulfur
        acyl_fragment = Chem.PathToSubmol(mol, acyl_atom_indices[2:])
        # Count carbon atoms in the fragment
        c_count = sum(1 for atom in acyl_fragment.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count < 4:
            return False, "Fatty acyl chain too short"
    else:
        return False, "Fatty acyl chain not found"

    # Check for phosphate groups (which may be deprotonated)
    phosphate_smarts = Chem.MolFromSmarts('OP(O)([O-])=O')
    num_phosphates = len(mol.GetSubstructMatches(phosphate_smarts))
    if num_phosphates < 3:
        return False, f"Found {num_phosphates} phosphate groups, expected at least 3"

    return True, "Molecule is a 3-hydroxy fatty acyl-CoA"

__metadata__ = {
    'chemical_class': {
        'name': '3-hydroxy fatty acyl-CoA(4-)',
        'definition': 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.',
    },
    'success': True,
    'error': '',
}
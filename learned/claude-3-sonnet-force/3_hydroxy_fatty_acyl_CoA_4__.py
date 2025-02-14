"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:36056 3-hydroxy fatty acyl-CoA(4-)
An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

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

    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(OP(=O)(OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12)O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Find the 3-hydroxy fatty acid chain
    chain_pattern = Chem.MolFromSmarts("[C;H3]([C;H2])[C@H](O)[C@@](C)(C)C(=O)C[C@@H]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) == 0:
        return False, "No 3-hydroxy fatty acid chain found"

    # Check fatty acid chain length
    chain_length = sum(1 for atom in mol.GetAtomWithIdx(chain_matches[0][3]).GetNeighbors() if atom.GetAtomicNum() == 6)
    if chain_length < 4 or chain_length > 24:
        return False, f"Fatty acid chain length ({chain_length}) not within expected range (4-24)"

    # Check for double bonds in the fatty acid chain
    double_bond_pattern = Chem.MolFromSmarts("[C@H]=C[C@H]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if double_bond_matches:
        # Ensure double bonds follow expected pattern
        allowed_double_bonds = [4, 7, 10, 13, 16, 19, 22]
        double_bond_positions = [match[1] for match in double_bond_matches]
        for pos in double_bond_positions:
            if pos not in allowed_double_bonds:
                return False, f"Double bond at position {pos} not allowed"

    # Calculate charge based on deprotonated phosphate and diphosphate groups
    deprotonated_groups = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1 and atom.GetAtomicNum() == 8)
    charge = -deprotonated_groups

    if charge != -4:
        return False, f"Expected charge of -4, found {charge}"

    return True, "Contains 3-hydroxy fatty acid chain with CoA group and expected charge"
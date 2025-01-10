"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester bond pattern: C(=O)S
    thioester_bond = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_bond):
        return False, "No thioester bond found"

    # Check for Coenzyme A structure: Contains adenosine triphosphate portion
    coA_structure = Chem.MolFromSmarts("NC(=O)C[C@H](O)C(C)(C)CO[P]([O-])(=O)O[P]([O-])(=O)OC...") # Truncated for simplicity
    if not mol.HasSubstructMatch(coA_structure):
        return False, "No Coenzyme A structure found"

    # Get the fatty acid chain attached to thioester bond (R in RC(=O)S)
    r_group_pattern = Chem.MolFromSmarts("S[C](=O)[C]-[CH2]-[CH2]-[CH2]")  # pattern indicates beginning of acyl chain
    for match in mol.GetSubstructMatches(r_group_pattern):
        chain_length = 0
        for idx in range(match[-1], mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon atom
                chain_length += 1
            else:
                break
        if 6 <= chain_length <= 12:
            break
    else:
        return False, f"Fatty acid chain not within medium-chain length, found length {chain_length}"

    # Check for deprotonated phosphate groups (count -P(=O)([O-])O patterns)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 2:
        return False, f"Insufficient deprotonated phosphate groups, found {phosphate_count}"

    # Ensure total charge is consistent with (4-)
    total_charge = rdMolDescriptors.CalcExactMolWt(mol, onlyHeavy=True, includeImplicitHs=False)
    if total_charge != -4:
        return False, f"Total charge is not -4"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"
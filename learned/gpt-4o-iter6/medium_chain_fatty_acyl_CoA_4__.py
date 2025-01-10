"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the thioester bond pattern: non-H atoms in C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if thioester_pattern is None or not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"

    # Check for presence of CoA structure using a truncated pattern
    coA_pattern = Chem.MolFromSmarts("NC(=O)C[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)O[C@H]1O[C@H](N)
                                      [C@@H](O[C@H]-*)[C@@H]1OP(=O)([O-])[O-]") # CoA key pattern
    if coA_pattern is None or not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A structure found"

    # Find the fatty acyl chain length, should be between 6 and 12 carbons
    chain_length = 0
    for match in mol.GetSubstructMatches(thioester_pattern):
        for idx in range(match[-1] + 1, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # carbon atom
                chain_length += 1
            else:
                break
    if not (6 <= chain_length <= 12):
        return False, f"Fatty acid chain not within medium-chain length, found length {chain_length}"

    # Ensure the molecule is deprotonated to a 4- charge (using phosphate pattern count)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))

    if phosphate_count < 3:
        return False, f"Insufficient deprotonated phosphate groups (expected 3-4), found {phosphate_count}"

    # Validate the formal charge of the molecule
    total_charge = Descriptors.CalcFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is not -4, found {total_charge}"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"
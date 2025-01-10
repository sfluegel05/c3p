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

    # Check for the thioester bond pattern C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"

    # Check for the CoA structure using an accurate pattern, capturing critical components
    coA_pattern = Chem.MolFromSmarts("CCC(=O)NCCSC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)O[C@H]1O[C@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H](C1)OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No complete Coenzyme A structure found"

    # Calculate the fatty acid chain length between the thioester carbonyl group and the rest of the chain
    matches = mol.GetSubstructMatches(thioester_pattern)
    for match in matches:
        chain_length = 0
        for idx in range(match[-1] + 1, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # Must be a carbon atom for chain confirmation
                break
            chain_length += 1
        if 6 <= chain_length <= 12:
            break
    else:
        return False, "Fatty acid chain not within medium-chain length"

    # Confirm the deprotonated state to carry a -4 charge
    total_charge = Descriptors.CalcFormalCharge(mol)
    if total_charge != -4:
        return False, f"Total charge is not -4, found {total_charge}"

    return True, "Molecule matches medium-chain fatty acyl-CoA(4-) structure"
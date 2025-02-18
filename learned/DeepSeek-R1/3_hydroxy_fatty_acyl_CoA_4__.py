"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    The structure must contain a 3-hydroxy acyl thioester group, three phosphate groups with four deprotonated oxygens,
    and have an overall charge of -4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for 3-hydroxy acyl thioester pattern (S-C(=O)-C-C-OH)
    thioester_pattern = Chem.MolFromSmarts("[S]-C(=O)-C-C(-O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Missing 3-hydroxy acyl thioester group"

    # Find all phosphorus atoms
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) != 3:
        return False, f"Found {len(p_atoms)} phosphorus atoms, need exactly 3"

    # Count deprotonated oxygens attached to phosphorus
    o_minus_count = 0
    for p in p_atoms:
        for neighbor in p.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                o_minus_count += 1
    if o_minus_count != 4:
        return False, f"Found {o_minus_count} deprotonated phosphate oxygens, need exactly 4"

    # Verify overall charge state
    if Chem.GetFormalCharge(mol) != -4:
        return False, f"Formal charge is {Chem.GetFormalCharge(mol)}, should be -4"

    return True, "Contains 3-hydroxy acyl thioester, three phosphates with four deprotonated oxygens, and charge -4"
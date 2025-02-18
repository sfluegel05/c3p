"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:XXXXX 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES.
    Must have: CoA structure, thioester linkage, and 3-hydroxy group on the fatty acid.

    Args:
        smiles (str): Input SMILES

    Returns:
        bool: True if matches criteria
        str: Reason for decision
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for CoA substructure (pantetheine-phosphate-ribose-adenine core)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing CoA substructure"

    # Find thioester linkage (S connected to carbonyl group)
    thioester = Chem.MolFromSmarts("[SX2][CX3](=O)")
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester bond"

    # Check 3-hydroxy position: -O-C-C-C(=O)S- (hydroxyl on third carbon from thioester)
    hydroxy_pattern = Chem.MolFromSmarts("[CX3](=O)-[SX2]-[CX4H2]-[CX4H2]-[CX3H1]([OX2H])")
    matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not matches:
        return False, "3-hydroxy group not found in correct position"

    # Verify the hydroxyl is on the beta carbon (third from thioester)
    # Additional check for chain length (at least 4 carbons including the thioester)
    chain_check = Chem.MolFromSmarts("[CX3](=O)-[SX2]-[CX4H2]-[CX4H2]-[CX3H1]([OX2H])-*")
    if not mol.HasSubstructMatch(chain_check):
        return False, "Insufficient chain length or hydroxyl position"

    return True, "3-hydroxy fatty acyl-CoA structure confirmed"
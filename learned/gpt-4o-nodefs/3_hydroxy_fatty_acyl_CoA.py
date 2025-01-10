"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    This is defined by the presence of a 3-hydroxy group on the fatty acid and the CoA structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 3-hydroxy group within a larger fatty acid chain 
    # -C(3)-[O]-C(=O)- can be more flexibly aware of longer carbon chains
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[CH2][CH2][CX3](=O)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No appropriate 3-hydroxy group in fatty acid chain found"
    
    # Coenzyme A moiety detection
    # Simplify or broaden pattern; CoA is typically a phosphate, pantetheine structure 
    coa_pattern = Chem.MolFromSmarts("C1[C@H](O)CO[C@H]1OP(=O)(O)OC[C@H]2[C@H](O)[C@H](O)[C@H](O2)OP(=O)(O)OCC[NH2+]CCCC(=O)[C@H](O)CCS")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    return True, "Contains 3-hydroxy fatty acid group with Coenzyme A moiety"
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
    
    # Check for 3-hydroxy fatty acid part: -C(3)-[O]-C(=O)-
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C[C@H](O)C(=O)")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid group found"
        
    # Check for CoA moiety: Coenzyme A structure recognition
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)OC[C@H]1O[C@@H]([C@@H](O)[C@H]1O)OP(=O)(O)OCSc1ncnc2n(ccn12)[C@H]3O[C@H]([C@H]([C@@H]3O)OP(=O)(O)O)CO")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    return True, "Contains 3-hydroxy fatty acid group with Coenzyme A moiety"
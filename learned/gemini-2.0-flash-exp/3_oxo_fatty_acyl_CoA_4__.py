"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
         str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core CoA pattern including diphosphate and terminal phosphate (deprotonated)
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([OX1,OX2])(=O)OP([OX1,OX2])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([OX1,OX2])([OX1,OX2])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
         return False, "Core CoA structure not found"
    
    # Check for the thioester linkage, and the adjacent carbon
    thioester_pattern = Chem.MolFromSmarts("[CX4](=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"
    
    #Check for 3-oxo group adjacent to a methylene group and the fatty acid chain
    oxo_group_pattern = Chem.MolFromSmarts("[CX4](=O)[CX2][CX4](=O)~[CX4]") # The 3-oxo group is next to a CH2 and also connects to a carbon of a chain
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "3-oxo group and adjacent methylene not found or improperly connected to fatty acid"

    # Check for a long aliphatic chain (at least 2 carbons) directly attached to the carbonyl
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4](=O)~[CX4][CX4]") # At least 2 carbons attached to the carbonyl of the 3-oxo group
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
      return False, "Missing or too short fatty acyl chain attached to 3-oxo group"


    # Check for at least 30 heavy atoms and 3 phosphorus atoms
    heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    if heavy_atoms < 30:
        return False, "Not enough heavy atoms."

    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 3:
        return False, "Not exactly 3 phosphorus atoms."
    
    o_minus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if o_minus_count < 4:
        return False, f"Must have at least 4 deprotonated oxygens."


    return True, "Meets all criteria for a 3-oxo-fatty acyl-CoA(4-)"
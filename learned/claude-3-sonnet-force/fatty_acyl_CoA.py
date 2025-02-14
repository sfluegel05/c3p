"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:35620 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA is an acyl-CoA that results from the formal condensation of the
    thiol group of coenzyme A with the carboxy group of any fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA substructure
    coa_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H]([C@H](OP(O)(O)=O)OP(O)(O)=O)O1)COP(O)(=O)OP(O)(=O)OC[C@H](C(=O)NCCC(=O)NCCS)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA substructure found"
    
    # Find ester functional group
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    for match in ester_matches:
        ester_carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in ester_carbon.GetNeighbors():
            if mol.HasSubstructMatch(fatty_acid_pattern, beginAtomIdx=neighbor.GetIdx()):
                break
        else:
            return False, "No fatty acid chain found"
    
    # Check molecular weight (typically >500 Da for fatty acyl-CoA)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for fatty acyl-CoA"
    
    return True, "Contains CoA substructure with a fatty acid chain attached via an ester linkage"
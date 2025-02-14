"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: CHEBI:36795 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is defined as a hydroxy fatty acyl-CoA that results from the
    formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

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

    # Check for coenzyme A moiety
    coa_pattern = Chem.MolFromSmarts("C(COP(OP(OC[C@H]1O[C@@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)(=O)O)(=O)O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing coenzyme A moiety"
    
    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4H2]([CX4H2])([CX4H2])([CX3](=[OX1])[CH3])[OX2H,OX1-]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"
    
    # Check for 3-hydroxy group on fatty acid chain
    for match in fatty_acid_matches:
        fatty_acid_chain = mol.GetAtomWithIdx(match[0])
        if fatty_acid_chain.GetProp("nHydroxyls") == "1":
            hydroxy_carbon = fatty_acid_chain.GetNeighbors()[0]
            if hydroxy_carbon.GetProp("nHydroxyls") == "1" and hydroxy_carbon.GetDegree() == 4:
                break
    else:
        return False, "No 3-hydroxy group found on fatty acid chain"

    # Count carbons and check molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Fatty acid chain too short"
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1500:
        return False, "Molecular weight outside expected range for 3-hydroxy fatty acyl-CoA"

    return True, "Contains a 3-hydroxy fatty acid chain connected to coenzyme A via a thioester bond"
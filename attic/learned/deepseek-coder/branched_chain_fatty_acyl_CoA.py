"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:28494 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA with a branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety using a simplified pattern
    # CoA has a characteristic adenosine diphosphate pantetheine structure
    coa_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][SX2][CX3](=O)[CX4]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Find the fatty acid part by looking for the carbonyl carbon attached to CoA
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2][CX4][CX4][NX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid part found"

    # Get the carbon chain attached to the carbonyl
    carbon_chain = []
    for match in fatty_acid_matches:
        carbonyl_carbon = match[0]
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in match:
                carbon_chain.append(neighbor.GetIdx())
                break

    # Check if the carbon chain is branched
    is_branched = False
    for atom_idx in carbon_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetDegree() > 2:  # Branch point
            is_branched = True
            break

    if not is_branched:
        return False, "Fatty acid chain is not branched"

    # Check molecular weight - should be >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for branched-chain fatty acyl-CoA"

    return True, "Contains CoA moiety with branched fatty acid chain"
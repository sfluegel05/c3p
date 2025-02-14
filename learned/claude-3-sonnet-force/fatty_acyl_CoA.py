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
    coa_pattern = Chem.MolFromSmarts("[C@H]1([NH2]C2=C(C(=NC=N2)N)N)O[C@@H]([C@@H]([C@@H]1OP(O)(O)=O)OP(O)(O)=O)COP(O)(=O)OP(O)(=O)OC[C@@H](C(=O)NCCC(=O)NCCS)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA substructure found"
    
    # Find ester carbonyl carbon
    ester_carbonyl = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetExplicitValence() == 3 and atom.GetHybridization() == Chem.HybridizationType.SP2:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == 0:
                    ester_carbonyl = atom
                    break
    
    if ester_carbonyl is None:
        return False, "No ester carbonyl found"
    
    # Recursively check for fatty acid chain
    def is_fatty_acid_chain(atom, min_length):
        chain_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetFormalCharge() == 0:
                chain_length += 1 + is_fatty_acid_chain(neighbor, min_length - 1)
        return chain_length >= min_length
    
    if not is_fatty_acid_chain(ester_carbonyl, 6):
        return False, "Fatty acid chain too short"
    
    return True, "Contains CoA substructure with a fatty acid chain attached via an ester linkage"

# Example usage
example_smiles = "CCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
result, reason = is_fatty_acyl_CoA(example_smiles)
print(f"Result: {result}, Reason: {reason}")
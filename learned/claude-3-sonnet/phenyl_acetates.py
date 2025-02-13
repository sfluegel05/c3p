"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
           with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a phenyl acetate and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of phenyl ring
    phenyl_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl ring found"

    # Check for acetate group attached to phenyl ring
    # [c] represents aromatic carbon, [OX2] is oxygen with 2 connections
    # [CX3](=[OX1])[CH3] represents the acetate group (C(=O)CH3)
    acetate_on_phenyl_pattern = Chem.MolFromSmarts("[c][OX2][CX3](=[OX1])[CH3]")
    matches = mol.GetSubstructMatches(acetate_on_phenyl_pattern)
    
    if not matches:
        return False, "No acetate group attached to phenyl ring"

    # Additional validation to ensure the oxygen is directly connected to the phenyl ring
    for match in matches:
        phenyl_carbon = mol.GetAtomWithIdx(match[0])
        if phenyl_carbon.IsInRing() and phenyl_carbon.GetIsAromatic():
            return True, "Contains phenyl ring with acetate group attached"
            
    return False, "Acetate group not properly connected to phenyl ring"

def test_examples():
    """Test function with known examples"""
    examples = [
        "CC(=O)Oc1ccccc1",  # phenyl acetate
        "COc1cc(COC(C)=O)ccc1OC(C)=O",  # vanillyl alcohol diacetate
        "COC(=O)c1ccc(OC(C)=O)cc1",  # methyl-4-acetoxybenzoate
        "CC(=O)Oc1cccc(O)c1",  # resorcinol monoacetate
        "CC(=O)Oc1ccc(C=O)cc1",  # 4-acetoxy benzaldehyde
        # Counter-examples
        "CC(=O)OCCC",  # non-phenyl acetate
        "c1ccccc1O",  # phenol without acetate
    ]
    
    for smi in examples:
        result, reason = is_phenyl_acetates(smi)
        print(f"SMILES: {smi}")
        print(f"Is phenyl acetate: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()
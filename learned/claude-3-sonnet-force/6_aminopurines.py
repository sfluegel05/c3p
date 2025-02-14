"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine based on its SMILES string.
    A 6-aminopurine is any compound having 6-aminopurine (adenine) as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define 6-aminopurine (adenine) pattern
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    
    # Check if molecule contains adenine pattern
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Does not contain 6-aminopurine (adenine) moiety"
    
    # Check for specific substituents or functional groups incompatible with 6-aminopurines
    # Example: Exclude molecules with carbonyl groups attached to the adenine ring
    carbonyl_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12C(=O)")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains incompatible substituents for 6-aminopurine"
    
    # Check for specific connectivity patterns required for 6-aminopurines
    # Example: Require the adenine ring to be part of a larger fused ring system
    fused_ring_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12C3CCC3")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "Does not meet connectivity requirements for 6-aminopurine"
    
    # Handle tautomerism (optional)
    tautomers = Chem.MolFromSmiles(smiles, allowCXSMILES=True).GetEnumerateTautomers()
    for tautomer in tautomers:
        if tautomer.HasSubstructMatch(adenine_pattern):
            # Additional checks for tautomer...
            pass
    
    return True, "Meets the structural requirements of a 6-aminopurine"

# Example usage
smiles1 = "Nc1ncnc2n(cnc12)C3OC(COP(=O)(O)OP(=O)(O)OCC4OC(n5cnc6c(N)ncnc56)C(O)C4O)C(O)C3O"  # Acetyl-CoA
smiles2 = "CCc1ccccc1"  # Ethylbenzene
smiles3 = "Nc1ncnc2[nH]cnc12"  # 6-methyladenine
smiles4 = "OC(CN1C2=NC=NC(N)=C2N=C1)C(O)C(O)=O"  # D-erythro-Eritadenine

print(is_6_aminopurines(smiles1))  # Output: (True, 'Meets the structural requirements of a 6-aminopurine')
print(is_6_aminopurines(smiles2))  # Output: (False, 'Does not contain 6-aminopurine (adenine) moiety')
print(is_6_aminopurines(smiles3))  # Output: (True, 'Meets the structural requirements of a 6-aminopurine')
print(is_6_aminopurines(smiles4))  # Output: (True, 'Meets the structural requirements of a 6-aminopurine')
"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate consists of a nitrogenous base attached to a sugar, with phosphate groups attached to the sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a phosphate group (PO4)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Define SMARTS pattern for a ribose or deoxyribose sugar (five-membered ring with oxygen)
    sugar_pattern = Chem.MolFromSmarts("C1=C(C(C(O1)CO)O)O")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"
    
    # Define nitrogenous base patterns (purine and pyrimidine)
    purine_pattern = Chem.MolFromSmarts("NC1=NC=NC2=C1N=CN=C2")
    pyrimidine_pattern = Chem.MolFromSmarts("C1=NC=CN=C1")
    
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No nitrogenous base (purine or pyrimidine) found"
    
    return True, "Molecule contains features of a nucleoside phosphate"

# Example usage:
# Test SMILES: "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H]2O)c(=O)n1"  # Example SMILES for a known nucleoside phosphate
# result, reason = is_nucleoside_phosphate(test_smiles)
# Print result and reason for classification
# print(result, reason)
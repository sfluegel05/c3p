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
        bool: True if the molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a phosphate group
    # Considering both organic phosphates (P(=O)(O)OP) and variations like cyclic, linked forms
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Define a more inclusive SMARTS pattern for a sugar (five-membered ring with an oxygen member)
    sugar_pattern = Chem.MolFromSmarts("C1OC(CO)C1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar pattern found"
    
    # Broaden nitrogenous base patterns to include more diverse structures
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]cnc12") 
    pyrimidine_pattern = Chem.MolFromSmarts("c1[nH]cnc1")
    
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No nitrogenous base (purine or pyrimidine) found"
    
    return True, "Molecule contains features of a nucleoside phosphate"

# Example usage:
# test_smiles = "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](OP(O)(O)=O)[C@H]2O)c(=O)n1" # Example SMILES for a known nucleoside phosphate
# result, reason = is_nucleoside_phosphate(test_smiles)
# print(result, reason)
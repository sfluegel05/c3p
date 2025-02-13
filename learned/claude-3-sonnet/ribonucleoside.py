"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:18021 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside is a nucleoside where the sugar component is D-ribose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for D-ribose sugar pattern
    ribose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No D-ribose sugar found"
    
    # Look for nucleobases (pyrimidine or purine)
    pyrimidine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(=O)[nH]c1")
    purine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    if not (mol.HasSubstructMatch(pyrimidine_pattern) or mol.HasSubstructMatch(purine_pattern)):
        return False, "No nucleobase found"
    
    # Look for N-glycosidic bond between ribose and nucleobase
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O)[nr3r5,nr5r5,nr5r6]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No N-glycosidic bond found"
    
    # Check for common functional groups and substituents
    has_hydroxyl = any(atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    has_phosphate = any(atom.GetSymbol() == 'P' for atom in mol.GetAtoms())
    
    # Check molecular weight and atom counts
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    num_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    num_n = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N')
    num_o = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    if not (200 < mol_wt < 600):
        return False, "Molecular weight out of typical range for ribonucleosides"
    if not (5 <= num_c <= 20):
        return False, "Carbon atom count out of typical range for ribonucleosides"
    if not (1 <= num_n <= 8):
        return False, "Nitrogen atom count out of typical range for ribonucleosides"
    if not (3 <= num_o <= 10):
        return False, "Oxygen atom count out of typical range for ribonucleosides"
    
    # If all checks pass, classify as a ribonucleoside
    reason = "Contains D-ribose sugar, nucleobase, and N-glycosidic bond"
    if has_hydroxyl:
        reason += ", with hydroxyl group(s)"
    if has_phosphate:
        reason += ", with phosphate group(s)"
    
    return True, reason
"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: CHEBI:33115 butyrate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is a carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester functional group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check if the ester group is part of a butyrate
    butyrate_pattern = Chem.MolFromSmarts("CCCC(=O)O[C]")
    butyrate_matches = mol.GetSubstructMatches(butyrate_pattern)
    if not butyrate_matches:
        return False, "Ester group is not part of a butyrate"
    
    # Check for isotopic labeling or other modifications
    butyl_patterns = [Chem.MolFromSmarts(f"[{isotope}CCCC](=O)O[C]") for isotope in ["", "13C", "14C", "2H"]]
    for pattern in butyl_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a butyrate ester group (with potential isotopic labeling)"
    
    # Additional checks based on molecular properties
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for butyrate esters"
    
    n_rotatable = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too few rotatable bonds for a butyrate ester"
    
    return True, "Contains a butyrate ester group"
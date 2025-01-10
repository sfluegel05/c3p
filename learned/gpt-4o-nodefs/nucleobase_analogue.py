"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for purine or pyrimidine rings
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]ncnc12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncnc[nH]1")
    
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No purine or pyrimidine core found"

    # Common functional group patterns in nucleobase analogues
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    if mol.HasSubstructMatch(oxo_pattern):
        return True, "Contains common oxo functional group for nucleobase analogues"

    amino_pattern = Chem.MolFromSmarts("NC")
    if mol.HasSubstructMatch(amino_pattern):
        return True, "Contains common amino functional group for nucleobase analogues"

    # Check for halogen substitutions
    halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogen_pattern):
        return True, "Contains halogen substitution, common in nucleobase analogues"

    # If none of the modifications above are found
    return False, "Lacks modifications typical of nucleobase analogues"
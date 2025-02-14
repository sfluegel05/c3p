"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide is a cerebroside with a glucose head group and a sphingosine linked to a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced pattern for beta-D-glucosyl, accounting for different linkages
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H](O)[C@@H](O)[C@H]1O)CO")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucosyl unit found"
    
    # Refined pattern for sphingosine backbone including amide bond
    sphingosine_pattern = Chem.MolFromSmarts("NC(=O)C[C@@H](O)CO")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine structure with amide linkage found"
    
    # Improved detection for long aliphatic chains
    long_chain_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCCCC"),
        Chem.MolFromSmarts("C=C")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns):
        return False, "No suitable aliphatic chain (indicative of fatty acid) found"
    
    return True, "Contains beta-D-glucosyl unit, sphingosine backbone with amide linkage, and long fatty acid chain"
"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Updated pattern for beta-D-glucosyl unit with flexible stereochemistry
    glucose_pattern = Chem.MolFromSmarts("O[C@H]([C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)CO")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucosyl unit found"

    # Refined sphingosine backbone detection
    sphingosine_pattern = Chem.MolFromSmarts("NC(=O)C[C@@H](O)C=C or NC(=O)C[C@@H](O)CC")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine structure with amide linkage found"
    
    # Patterns to identify long aliphatic chains (both saturated and unsaturated)
    long_chain_saturated = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC")
    long_chain_unsaturated = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(long_chain_saturated) and not mol.HasSubstructMatch(long_chain_unsaturated):
        return False, "Suitable aliphatic chain (indicative of fatty acid) not satisfactorily found"
    
    return True, "Contains beta-D-glucosyl unit, sphingosine backbone with amide linkage, and long fatty acid chain"
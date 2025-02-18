"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants 
    and may vary in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Steroid backbone pattern: Tetracyclic
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1(CCC3C2CCC4C3(CCCC4)C)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No tetracyclic steroid backbone found"
    
    # Check for various common side chains and modifications
    # We'll use basic checks here - in practice phytosterols can vary widely
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(double_bond_pattern):
        shared_double_bonds = True
    else:
        shared_double_bonds = False
        
    # Represent common side chain variations for phytosterols as SMARTS patterns
    farnesyl_chain_pattern = Chem.MolFromSmarts("C(C)CC=C(C)C")
    if mol.HasSubstructMatch(farnesyl_chain_pattern):
        side_chain_match = True
    else:
        side_chain_match = False

    # Conclusion based on substructure matches
    if shared_double_bonds and side_chain_match:
        return True, "Contains common phytosterol modifications (double bonds and possible farnesyl chain)"
    
    return False, "Lacks distinct side chain modifications typical of phytosterols"
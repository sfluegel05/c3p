"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with acyl group at 1-position
    glycerol_acyl_pattern = Chem.MolFromSmarts("C(CO[P](=O)(O)OCCN)(OC(=O)[R])O") 
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "No glycerol-acyl pattern found at 1-position"

    # Look for phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("OP(=O)(OCCN)O")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "Phosphoethanolamine group not found"
    
    # Check for presence of long acyl chains via carbon count, associated with fatty acids
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 15:
        return False, f"Not enough carbon atoms for a typical acyl chain, found {num_carbons}"

    return True, "Contains 1-O-acylglycerophosphoethanolamine structure with appropriate groups"
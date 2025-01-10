"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids include sphinganine, its homologs and stereoisomers, and hydroxyl or unsaturated derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for characteristic amino alcohol moiety in sphingoids
    amino_alcohol_pattern = Chem.MolFromSmarts("C[C@H](N)CO")  # More specific pattern for amino alcohol
    if not mol.HasSubstructMatch(amino_alcohol_pattern):
        return False, "No characteristic amino alcohol moiety found"
    
    # Check for a long hydrocarbon chain
    chain_pattern = Chem.MolFromSmarts("C" * 10)  # Minimum length of 10 carbon atoms
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Capture potential unsaturation (double bonds) in the hydrocarbon chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No unsaturation found, but this is optional"

    # Check for stereochemistry in the hydrocarbon portion
    stereo_present = any(atom.HasProp('_CIPCode') for atom in mol.GetAtoms())  # Identifies stereocenters
    if not stereo_present:
        return False, "No stereochemistry typical of sphingoids found"

    return True, "SMILES corresponds to a recognized sphingoid structure"
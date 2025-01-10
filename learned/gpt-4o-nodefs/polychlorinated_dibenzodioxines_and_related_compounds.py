"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule belongs to the class of polychlorinated dibenzodioxines and related compounds
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class description, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for chlorine atoms
    cl_count = sum(atom.GetSymbol() == 'Cl' for atom in mol.GetAtoms())
    if cl_count < 2:
        return False, "Less than two chlorine atoms"

    # Check for biphenyl structure (two phenyl rings connected)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if mol.HasSubstructMatch(biphenyl_pattern):
        return True, "Biphenyl structure with chlorines found"

    # Check for dibenzodioxin structure: two benzene rings with two oxygens
    dibenzodioxin_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3oc2cc1")
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        return True, "Dibenzodioxin structure with chlorines found"
    
    # Check for dibenzofuran structure
    dibenzofuran_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3c2cc1")
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        return True, "Dibenzofuran structure with chlorines found"

    return False, "Does not match any known structure for this class"
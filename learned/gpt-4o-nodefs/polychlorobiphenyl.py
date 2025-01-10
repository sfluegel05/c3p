"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl (PCB) based on its SMILES string.
    PCBs have a biphenyl backbone with multiple chlorine substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PCB, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for biphenyl pattern: two connected benzene rings
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl backbone found"
        
    # Count the number of chlorine atoms
    cl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
        
    if cl_count > 2:
        return True, f"Contains biphenyl backbone with {cl_count} chlorine substitutions"
    else:
        return False, f"Contains biphenyl backbone but only {cl_count} chlorine substitutions, need more than 2"

__metadata__ = {'chemical_class': {'id': 'None',
                                   'name': 'polychlorobiphenyl',
                                   'definition': 'Polychlorinated biphenyls (PCBs) are a group '
                                                 'of man-made organic chemicals consisting of carbon, '
                                                 'hydrogen and chlorine atoms. PCBs have a biphenyl '
                                                 'backbone with chlorine substitutions.'}}
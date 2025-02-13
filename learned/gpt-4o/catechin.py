"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Aiming to capture the flavan-3-ol skeleton with additional common substitutions,
    which includes two phenolic rings connected by a pyran ring with hydroxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broad SMARTS pattern for a flavan-3-ol core with possible variations
    flavan3ol_pattern = Chem.MolFromSmarts("C1=C(C=C(C=C1)O)C2(O)C3=C(C=C(C=C3)O)OC(CC2)C")
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol skeleton found"
    
    # Verify presence of hydroxyl or methoxy like groups  
    hydroxyl_or_methoxy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if hydroxyl_or_methoxy < 4:  # Typically a catechin has at least 4 oxygens, either in -OH or -OCH3
        return False, f"Insufficient -OH/-OCH3 groups, found {hydroxyl_or_methoxy}"
    
    # Ensure there are chiral centers as catechins are chiral molecules
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not stereo_centers:
        return False, "No stereochemistry observed, unlikely to be a catechin"
    
    return True, "Molecule matches flavan-3-ol skeleton with common hydroxyl/methoxy groups"
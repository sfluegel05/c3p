"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins typically have a flavan-3-ol skeleton which consists of two phenolic rings
    connected by a chiral three-carbon unit forming a pyran ring, often with hydroxyl and
    methyl ether groups.
    
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
    
    # Define a SMARTS pattern for a core flavan-3-ol. Consider some substitution variations.
    flavan3ol_pattern = Chem.MolFromSmarts("c1cc(O)ccc1-C2C[C@H](O)Oc3c(O)cc(O)cc23")
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol skeleton found"
    
    # Check for multiple hydroxyl groups beyond the core requirement
    # Catechins often have additional -OH groups or their equivalents (e.g., methoxy).
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    if hydroxyl_count < 4:  # Adjust based on expected range
        return False, f"Too few hydroxyl or equivalent groups, found {hydroxyl_count}"
    
    # Check for stereo configurations typically present in catechins
    stereo_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not stereo_centers:
        return False, "No chiral centers found for a catechin structure"
    
    return True, "Molecule matches flavan-3-ol skeleton with appropriate substitutions"
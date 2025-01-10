"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol typically contains a characteristic dihydroflavanone core structure
    with hydroxylation and specific chiral centers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjust the SMARTS pattern to detect a broader range of dihydroflavonol structures
    flavanone_pattern = Chem.MolFromSmarts("[C@H]1(C=O)Oc2cc(O)cc(O)c2C1")  # Basic flavanone core
    if not mol.HasSubstructMatch(flananone_pattern):
        return False, "No dihydroflavanone core found"
    
    # Check for necessary chiral centers, potential fine-tuning required for stereo-specific compounds
    stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(stereo) < 2:
        return False, "Chiral centers do not match typical dihydroflavonol structure"
    
    # Look for hydroxylation on the phenyl group and the pyran ring - basic pattern coverage
    hydroxyl_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Dihydroflavonol likely incomplete, lacking expected hydroxylation"
    
    # Validate that it matches known hydroxylated dihydroflavonol structures
    n_hydroxyls = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetSymbol() == 'O')
    if n_hydroxyls < 3:  # This threshold can be adjusted if needed
        return False, f"Found {n_hydroxyls} hydroxyl groups, expected at least 3"
    
    return True, "Contains dihydroflavonol features including core structure, chiral centers, and appropriate hydroxylation"
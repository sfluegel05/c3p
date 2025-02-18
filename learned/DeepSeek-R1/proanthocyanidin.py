"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: CHEBI:26154 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    Proanthocyanidins are flavonoid oligomers formed by condensation of hydroxyflavan units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavan unit check - look for benzopyran structure with hydroxyls
    # Flavan pattern: benzene ring (A) connected to dihydropyran (C) with oxygen
    flavan_core = Chem.MolFromSmarts("[c]1[c][c][c]([OH])[c][c]1-C1-C(-O-)C-C-C1")
    if not mol.HasSubstructMatch(flavan_core):
        return False, "No flavan core structure found"

    # Find all flavan units (benzopyran with at least two hydroxyls)
    # A-ring: at least two hydroxyls; C-ring: pyran oxygen
    flavan_pattern = Chem.MolFromSmarts(
        "[OH]c1c([OH])ccc2C(C)(O)CCCC12")
    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    
    if len(flavan_matches) < 2:
        return False, f"Found {len(flavan_matches)} flavan units, need at least 2"

    # Check for interflavan linkages (direct C-C bonds between aromatic carbons)
    # Look for two aromatic carbons connected by a single bond (C4-C8 type)
    linkage_pattern = Chem.MolFromSmarts("[c]-[CX4;H0]-[c]")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No interflavan C-C linkages detected"

    # Verify oligomer characteristics through molecular weight (dimers start ~500)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight {mol_wt:.1f} too low for oligomer"

    # Check for presence of galloyl groups (common but not mandatory)
    galloyl = Chem.MolFromSmarts("c1c(O)c(O)c(O)cc1C(=O)O")
    if galloyl and mol.HasSubstructMatch(galloyl):
        return True, "Contains flavan oligomer with galloyl group"

    return True, "Contains multiple flavan units with C-C interflavan linkages"
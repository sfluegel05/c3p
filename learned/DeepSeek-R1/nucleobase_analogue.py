"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogues
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue has a core structure similar to natural nucleobases (adenine, guanine, cytosine, thymine/uracil)
    but with chemical modifications.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define SMARTS patterns for natural nucleobase cores
    nucleobase_patterns = {
        "adenine": Chem.MolFromSmarts("Nc1ncnc2ncnc12"),
        "guanine": Chem.MolFromSmarts("Nc1nc2c(ncn2)c(=O)[nH]1"),
        "cytosine": Chem.MolFromSmarts("Nc1ncc(=O)[nH]c1=O"),
        "thymine/uracil": Chem.MolFromSmarts("O=C1NC(=O)CC=C1")
    }
    
    # Check if the molecule matches any nucleobase core
    core_matched = None
    for name, pattern in nucleobase_patterns.items():
        if mol.HasSubstructMatch(pattern):
            core_matched = name
            break
    
    if not core_matched:
        return False, "No nucleobase core detected"
    
    # Check for modifications: presence of substituents beyond the core
    # Get atoms in the core pattern
    core_atoms = set()
    for match in mol.GetSubstructMatches(pattern):
        core_atoms.update(match)
    
    # Check for non-core atoms connected to the core
    has_modification = False
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        # If one atom is in core and the other isn't, it's a substituent
        if (begin_idx in core_atoms) != (end_idx in core_atoms):
            has_modification = True
            break
    
    if not has_modification:
        return False, f"Matched {core_matched} core but no modifications"
    
    # Exclude exact natural nucleobases
    natural_smiles = {
        "adenine": "Nc1ncnc2ncnc12",
        "guanine": "Nc1nc2c(ncn2)c(=O)[nH]1",
        "cytosine": "Nc1cc(=O)[nH]c(=O)n1",
        "uracil": "O=C1NC(=O)CC=C1",
        "thymine": "Cc1ncc(=O)[nH]c1=O"
    }
    canon_smiles = Chem.CanonSmiles(smiles)
    if any(Chem.CanonSmiles(nat) == canon_smiles for nat in natural_smiles.values()):
        return False, "Exact natural nucleobase"
    
    return True, f"Modified {core_matched} core detected"
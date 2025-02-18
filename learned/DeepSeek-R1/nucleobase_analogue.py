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
    
    # Define SMARTS patterns for natural nucleobase cores (more general)
    nucleobase_patterns = {
        "adenine": Chem.MolFromSmarts("[NH2]C1=NC=NC2=C1N=CN2"),
        "guanine": Chem.MolFromSmarts("Nc1nc2c(ncn2)c(=O)[nH]1"),
        "cytosine": Chem.MolFromSmarts("Nc1nc(=O)ccn1"),
        "uracil": Chem.MolFromSmarts("O=C1NC(=O)C=C1"),
        "thymine": Chem.MolFromSmarts("CC1=CNC(=O)NC1=O")
    }
    
    # Check if the molecule matches any nucleobase core
    core_matched = None
    matched_pattern = None
    for name, pattern in nucleobase_patterns.items():
        if pattern and mol.HasSubstructMatch(pattern):
            core_matched = name
            matched_pattern = pattern
            break
    
    if not core_matched:
        return False, "No nucleobase core detected"
    
    # Get atoms in the matched core pattern
    core_atoms = set()
    for match in mol.GetSubstructMatches(matched_pattern):
        core_atoms.update(match)
    
    # Check for modifications: atoms outside core or different functional groups
    has_modification = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in core_atoms:
            has_modification = True
            break
        else:
            # Check if atom properties differ from core (e.g., different element, charge)
            # This part is complex; for simplicity, assume any substitution is a modification
            pass
    
    if not has_modification:
        return False, f"Matched {core_matched} core but no modifications"
    
    # Exclude exact natural nucleobases by canonical SMILES comparison
    natural_smiles = {
        "adenine": "Nc1ncnc2ncnc12",
        "guanine": "Nc1nc2c(ncn2)c(=O)[nH]1",
        "cytosine": "Nc1nc(=O)ccn1",
        "uracil": "O=C1NC(=O)C=C1",
        "thymine": "Cc1ncc(=O)[nH]c1=O"
    }
    
    # Precompute canonical SMILES of natural bases
    natural_canon = set()
    for name, smi in natural_smiles.items():
        mol_nat = Chem.MolFromSmiles(smi)
        if mol_nat:
            natural_canon.add(Chem.CanonSmiles(Chem.MolToSmiles(mol_nat)))
    
    # Canonicalize input
    input_canon = Chem.CanonSmiles(smiles)
    if input_canon in natural_canon:
        return False, "Exact natural nucleobase"
    
    return True, f"Modified {core_matched} core detected"
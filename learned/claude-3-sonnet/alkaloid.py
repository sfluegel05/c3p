"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:38145 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is defined as a naturally occurring, basic nitrogen compound (mostly heterocyclic)
    occurring mostly in the plant kingdom, but also found in bacteria, fungi, and animals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for basic nitrogen
    has_basic_nitrogen = any(atom.GetTotalNumHs() < 3 and atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_basic_nitrogen:
        return False, "No basic nitrogen present"
    
    # Check for heterocyclic rings
    ring_info = mol.GetRingInfo()
    has_heterocyclic_ring = any(ring.IsHeterocyclic() for ring in ring_info.AtomRings())
    if not has_heterocyclic_ring:
        return False, "No heterocyclic rings found"
    
    # Check for common alkaloid substructures
    alkaloid_patterns = [
        Chem.MolFromSmarts("[n;H0,H1]"),       # Pyridine-like N
        Chem.MolFromSmarts("[n;H0](:c):c"),    # Pyrrole-like N
        Chem.MolFromSmarts("[n;H1](:c)(:c)"),  # Imidazole-like N
        Chem.MolFromSmarts("[N+]"),            # Quaternary N
        Chem.MolFromSmarts("[N+]#[C-]"),       # Isonitrile
        Chem.MolFromSmarts("C1=CN=C2NCCCC12"), # Piperidine
        Chem.MolFromSmarts("C1=CN=CC2=C1NC=C2"),  # Quinoline
        Chem.MolFromSmarts("C1=CN=CN2CCCC=12")    # Imidazole
    ]
    has_alkaloid_substructure = any(mol.HasSubstructMatch(pattern) for pattern in alkaloid_patterns)
    
    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if mol_wt < 150 or mol_wt > 1000:
        return False, "Molecular weight outside typical alkaloid range (150-1000 Da)"
    
    if n_rotatable < 2:
        return False, "Too few rotatable bonds for an alkaloid"
    
    if has_alkaloid_substructure:
        return True, "Contains basic nitrogen, heterocyclic rings, and common alkaloid substructures"
    else:
        return True, "Contains basic nitrogen and heterocyclic rings, but no common alkaloid substructures"
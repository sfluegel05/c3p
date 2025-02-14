"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as any steroid lactone that is a C28 steroid 
    with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get Murcko scaffold (core steroid structure)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # Check number of rings in the scaffold
    ring_info = scaffold.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, "No steroid backbone found (requires at least 4 fused rings)"
    
    # Count total number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Too few carbons for a C28 steroid (found {c_count} carbons)"
    
    # Look for lactone ring (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C1=CC(=O)OC1")  # General lactone pattern
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if len(lactone_matches) == 0:
        return False, "No lactone ring found"

    # Ensure lactone ring is part of a side chain (not fused to steroid core)
    # Extract atoms from steroid scaffold
    scaffold_atoms = set(atom.GetIdx() for atom in scaffold.GetAtoms())
    lactone_in_side_chain = False
    for match in lactone_matches:
        # Check if lactone ring atoms are outside the steroid core
        if not scaffold_atoms.intersection(match):
            lactone_in_side_chain = True
            break
    if not lactone_in_side_chain:
        return False, "Lactone ring is not in the side chain"

    return True, "Molecule is a withanolide (C28 steroid with steroid backbone and side-chain lactone ring)"
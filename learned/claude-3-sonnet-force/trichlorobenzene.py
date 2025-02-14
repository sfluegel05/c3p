"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: CHEBI:32221 trichlorobenzene

Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if chlorine_count != 3:
        return False, f"Molecule contains {chlorine_count} chlorine atoms, should be 3"
    
    # Find aromatic rings
    aromatic_rings = Chem.GetAromaticRings(mol)
    
    # Check if there is an aromatic benzene ring with 3 chlorine substituents
    benzene_rings_with_3_chlorines = []
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        ring_chlorines = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17 and atom.GetIdx() in ring_atoms]
        if len(ring_atoms) == 6 and len(ring_chlorines) == 3:
            benzene_rings_with_3_chlorines.append(ring_atoms)
    
    if not benzene_rings_with_3_chlorines:
        return False, "No benzene ring with 3 chlorine substituents found"
    
    return True, "Molecule contains a benzene ring with three chlorine substituents"
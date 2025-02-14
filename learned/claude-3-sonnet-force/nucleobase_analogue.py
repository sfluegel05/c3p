"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:37585 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.AromaticRings import AromaticRings

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aromatic rings
    aromatic_rings = AromaticRings.IdentifyAromaticRings(mol)
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    # Check for planarity
    if not AllChem.PlanarityAtomsSelfIter(mol):
        return False, "Molecule is not planar"
    
    # Check for typical heteroatoms and functional groups
    allowed_atoms = set([6, 7, 8, 9, 16, 17])  # C, N, O, F, S, Cl
    atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    if not all(atom in allowed_atoms for atom in atoms):
        return False, "Molecule contains unusual atoms for a nucleobase analogue"
    
    # Check for typical size range
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    if n_heavy_atoms < 5 or n_heavy_atoms > 25:
        return False, "Size deviates too much from typical nucleobases"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 50 or mol_wt > 300:
        return False, "Molecular weight deviates too much from typical nucleobase analogues"
    
    # If all conditions are met, classify as nucleobase analogue
    return True, "Contains aromatic rings, is planar, and has typical structure and size of a nucleobase analogue"
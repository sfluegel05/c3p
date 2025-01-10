"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring based on its SMILES string.
    An azole is a five-membered aromatic heterocycle containing at least one nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an azole ring, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate 2D coordinates for ring detection
    AllChem.Compute2DCoords(mol)

    # Find all rings in the molecule
    rings = mol.GetRingInfo()
    
    # Check each ring of size 5
    for ring_atoms in rings.AtomRings():
        if len(ring_atoms) != 5:
            continue
            
        # Get the atoms in the ring
        ring_atom_objects = [mol.GetAtomWithIdx(i) for i in ring_atoms]
        
        # Check if ring contains at least one nitrogen
        has_nitrogen = False
        other_heteroatoms = False
        all_aromatic = True
        
        for atom in ring_atom_objects:
            # Check for nitrogen
            if atom.GetAtomicNum() == 7:
                has_nitrogen = True
            # Check for other heteroatoms
            elif atom.GetAtomicNum() in [8, 16]:  # O or S
                other_heteroatoms = True
            # Check if all atoms are aromatic
            if not atom.GetIsAromatic():
                all_aromatic = False
                
        # If we found a five-membered aromatic ring with at least one nitrogen
        if has_nitrogen and all_aromatic:
            msg = "Contains a five-membered aromatic ring with nitrogen"
            if other_heteroatoms:
                msg += " and additional heteroatoms"
            return True, msg

    return False, "No azole ring found"
"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemistry from the input molecule
    mol_no_stereo = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol_no_stereo)

    # Define the hopane molecule
    hopane_smiles = "CC1(C)CCC2(C)CCC3CC4(C)CCC5CC(C)(CCCC5(C)C)CC4(C)C3CC2C1"
    hopane_mol = Chem.MolFromSmiles(hopane_smiles)
    if hopane_mol is None:
        return False, "Failed to generate hopane structure"

    # Get ring atoms of hopane_mol
    ring_info = hopane_mol.GetRingInfo()
    ring_atom_idxs = set()
    for ring in ring_info.AtomRings():
        ring_atom_idxs.update(ring)

    # Create submol of hopane ring skeleton
    hopane_ring_mol = Chem.PathToSubmol(hopane_mol, sorted(ring_atom_idxs))
    
    # Remove hydrogens from ring mol
    hopane_ring_mol = Chem.RemoveHs(hopane_ring_mol)
    
    # Prepare input molecule for matching
    mol_no_stereo = Chem.RemoveHs(mol_no_stereo)
    
    # Perform substructure match to check for the hopane ring system
    if mol_no_stereo.HasSubstructMatch(hopane_ring_mol):
        return True, "Contains hopane ring skeleton characteristic of hopanoids"
    else:
        return False, "Does not contain hopane ring skeleton"
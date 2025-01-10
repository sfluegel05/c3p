"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a macrocyclic lactone
    ring_info = mol.GetRingInfo()
    found_macrocyclic_lactone = False
    for ring in ring_info.BondRings():
        # Check if the ring is large (12 or more atoms) and contains an ester linkage (O=C-O)
        if len(ring) >= 12:
            ring_atoms = [mol.GetAtomWithIdx(bond.GetBeginAtomIdx()) for bond in [mol.GetBondWithIdx(idx) for idx in ring]]
            # Identify ester linkages in the ring
            ester_pattern = Chem.MolFromSmarts('[C]=[O][O]')
            if any(atom.GetSmarts() == 'C(=O)O' for atom in ring_atoms):
                found_macrocyclic_lactone = True
                break
    
    if not found_macrocyclic_lactone:
        return False, "No macrocyclic lactone with 12 or more members found"

    # We could also further analyze here for polyketide patterns, but for simplicity, we assume any valid large lactone is a possible polyketide derivative

    return True, "Contains a macrocyclic lactone with 12 or more members"
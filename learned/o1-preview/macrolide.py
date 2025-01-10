"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
"""

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Descriptors

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

    # Find all rings in the molecule
    ssr = Chem.GetSymmSSSR(mol)
    macrocycle_found = False

    for ring in ssr:
        ring_size = len(ring)
        # Check for rings with 12 or more atoms
        if ring_size >= 12:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_bonds = []
            for atom in ring_atoms:
                for bond in atom.GetBonds():
                    if bond.GetBeginAtom() in ring_atoms and bond.GetEndAtom() in ring_atoms:
                        if bond not in ring_bonds:
                            ring_bonds.append(bond)
            # Create a molecule object of the ring
            ring_mol = Chem.PathToSubmol(mol, ring)
            # Look for lactone functional group within the ring (-C(=O)O-)
            lactone_pattern = Chem.MolFromSmarts("C(=O)O")
            if ring_mol.HasSubstructMatch(lactone_pattern):
                macrocycle_found = True
                break  # Exit after finding one macrocyclic lactone

    if not macrocycle_found:
        return False, "No macrocyclic lactone ring of 12 or more members found"

    return True, "Contains a macrocyclic lactone ring of 12 or more members"
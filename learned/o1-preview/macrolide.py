"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
"""

from rdkit import Chem

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

    # Iterate over each ring
    for ring in ssr:
        ring_size = len(ring)
        # Check for rings with 12 or more atoms
        if ring_size >= 12:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            ring_atom_idxs = set(ring)
            # Search for lactone functionality within the ring
            for atom in ring_atoms:
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check if carbon is a carbonyl carbon within the ring
                    is_carbonyl = False
                    for nbr in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 8:
                            is_carbonyl = True  # Found C=O
                            break
                    if is_carbonyl:
                        # Check for single bond to oxygen within the ring (ester linkage)
                        for nbr in atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nbr.GetAtomicNum() == 8:
                                # Ensure the oxygen is in the ring
                                if nbr.GetIdx() in ring_atom_idxs:
                                    macrocycle_found = True
                                    return True, "Contains a macrocyclic lactone ring of 12 or more members"
    if not macrocycle_found:
        return False, "No macrocyclic lactone ring of 12 or more members found"
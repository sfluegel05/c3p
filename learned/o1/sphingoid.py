"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: sphingoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is defined as sphinganine, its homologs and stereoisomers,
    and the hydroxy and unsaturated derivatives of these compounds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find primary amine nitrogen atoms (degree 2 nitrogen with no charge)
    amine_atoms = [atom for atom in mol.GetAtoms()
                   if atom.GetAtomicNum() == 7 and atom.GetDegree() == 2 and atom.GetFormalCharge() == 0]
    if not amine_atoms:
        return False, "No primary amino group found"

    for amine_atom in amine_atoms:
        # Get the carbon atom attached to the amine nitrogen (C2)
        neighbors = [nbr for nbr in amine_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not neighbors:
            continue  # No carbon attached to nitrogen
        c2_atom = neighbors[0]

        # Initialize variables for chain traversal
        visited = set()
        chain_atoms = []

        def traverse_chain(atom):
            if atom.GetAtomicNum() != 6:
                return
            if atom.GetIdx() in visited:
                return
            visited.add(atom.GetIdx())
            chain_atoms.append(atom)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).IsInRing() == False:
                    traverse_chain(nbr)

        # Traverse the chain starting from c2_atom
        traverse_chain(c2_atom)
        chain_length = len(chain_atoms)

        if chain_length < 12:
            continue  # Chain too short, try next amine

        # Check for hydroxyl groups attached to carbons adjacent to c2_atom
        has_hydroxyl = False
        for nbr in c2_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != amine_atom.GetIdx():
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 8 and nbr2.GetDegree() == 1:
                        has_hydroxyl = True
                        break
                if has_hydroxyl:
                    break
        if not has_hydroxyl:
            continue  # No hydroxyl group adjacent to amino group

        # If all checks passed
        return True, "Molecule matches sphingoid structural features"

    return False, "Molecule does not match sphingoid structural features"
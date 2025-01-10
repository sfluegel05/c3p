"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a flavanone with a hydroxy substituent at the 4' position of the B-ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core SMARTS pattern (corrected)
    flavanone_core_smarts = 'O=C1CCOc2ccccc12'  # Corrected flavanone core
    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if flavanone_core is None:
        return None, "Error in flavanone core SMARTS pattern"

    # Check if the molecule contains the flavanone core
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Flavanone core not found"

    # Define SMARTS pattern for hydroxy group at 4' position on B-ring (corrected)
    hydroxy_b_ring_smarts = 'O=C1CCOc2ccc(O)cc12'  # Corrected pattern with hydroxy on B-ring
    hydroxy_b_ring = Chem.MolFromSmarts(hydroxy_b_ring_smarts)
    if hydroxy_b_ring is None:
        return None, "Error in 4'-hydroxyflavanone SMARTS pattern"

    # Check if molecule has the hydroxy group at 4' position
    if mol.HasSubstructMatch(hydroxy_b_ring):
        return True, "Molecule is a 4'-hydroxyflavanone"
    else:
        # Further check for substituted B-rings with hydroxy at 4' position
        # Use a more general pattern for flavanone core with substitutions allowed
        flavanone_general_smarts = 'O=C1CCOc2ccccc12'  # General flavanone core
        flavanone_general = Chem.MolFromSmarts(flavanone_general_smarts)
        matches = mol.GetSubstructMatches(flavanone_general)
        if matches:
            for match in matches:
                # Get the atom indices of the B-ring carbons
                b_ring_atoms = [match[i] for i in range(6,12)]  # Atoms of B-ring in the match
                for atom_idx in b_ring_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    # Check for hydroxy group attached to B-ring carbons
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                            return True, "Molecule is a 4'-hydroxyflavanone with substituted B-ring"
            return False, "4'-hydroxy group on B-ring not found"
        else:
            return False, "Flavanone core not found"
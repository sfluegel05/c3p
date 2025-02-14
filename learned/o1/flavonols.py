"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: Flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone in which the hydrogen at position 3 of the heterocyclic ring
    is replaced by a hydroxy group (i.e., it's a 3-hydroxyflavone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core SMARTS pattern (benzopyran-4-one structure)
    flavone_smarts = 'c1cc2oc(=O)cc2c1'  # Flavone core structure
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    if flavone_pattern is None:
        return False, "Invalid SMARTS pattern for flavone core"

    # Find matches for the flavone core
    matches = mol.GetSubstructMatches(flavone_pattern)
    if not matches:
        return False, "Does not contain the flavone core structure"

    # For each match, check for hydroxyl at position 3 of the heterocyclic ring
    for match in matches:
        # Atom indices in the flavone pattern:
        # match[0] to match[6] corresponding to atoms in the SMARTS pattern

        # Position of the carbon at position 3 in the heterocyclic ring (pyran ring)
        # which is the 6th atom in the SMARTS pattern (0-based index)
        position3_atom_idx = match[5]
        position3_atom = mol.GetAtomWithIdx(position3_atom_idx)

        # Check if the position 3 carbon has a hydroxyl group attached
        has_OH = False
        for neighbor in position3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                bond = mol.GetBondBetweenAtoms(position3_atom_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if the oxygen is a hydroxyl group (i.e., has one hydrogen)
                    if neighbor.GetTotalNumHs() == 1:
                        has_OH = True
                        break
        if has_OH:
            return True, "Contains flavonol core with hydroxyl at position 3 characteristic of flavonols"
        else:
            continue

    return False, "Does not have hydroxyl group at position 3 of the flavone core"
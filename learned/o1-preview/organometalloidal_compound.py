"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound is defined as 'A compound having bonds between one or more metalloid atoms
    and one or more carbon atoms of an organyl group.'

    Metalloid elements considered here are arsenic (As) and antimony (Sb).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # List of metalloid atomic numbers to consider
    metalloid_atomic_nums = [33, 51]  # As, Sb

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to check for metalloid-carbon bond
    has_metalloid_carbon_bond = False

    # Iterate over metalloid atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloid_atomic_nums:
            metalloid_atom = atom
            # Check neighbors of metalloid atom
            for neighbor in metalloid_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    # Check if carbon is sp3-hybridized (non-aromatic, single bonds)
                    if (not neighbor.IsInRing() and 
                        neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and
                        not neighbor.GetIsAromatic()):
                        return True, f"Metalloid atom ({metalloid_atom.GetSymbol()}) bonded to sp3-hybridized carbon atom of organyl group"
                    else:
                        # Continue checking other neighbors
                        continue
    # If no suitable metalloid-carbon bonds found
    return False, "No metalloid-carbon bonds to sp3-hybridized carbon atoms of organyl groups found"
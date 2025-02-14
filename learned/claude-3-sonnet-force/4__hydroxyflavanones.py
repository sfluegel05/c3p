"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
"""
Classifies: CHEBI:28096 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a flavanone core with a hydroxy substituent located at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for flavanone core pattern (two fused rings with C=O in one ring)
    flavanone_pattern = Chem.MolFromSmarts("[C&r5,r6]1=C[C@H](=O)[C@]2=C[C@H]=C[C@H]=C12")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core found"

    # Get the atom at position 4' of the flavanone core
    flavor_core_ring = mol.GetRingInfo().BondRings()[0]
    pos4_atom_idx = flavor_core_ring[2]
    pos4_atom = mol.GetAtomWithIdx(pos4_atom_idx)

    # Check if the atom at position 4' has a hydroxy substituent
    if pos4_atom.GetTotalNumHs() == 1 and pos4_atom.GetAtomicNum() == 6:
        for neighbor in pos4_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                return True, "Contains flavanone core with a hydroxy group at position 4'"

    return False, "No hydroxy group found at position 4' of the flavanone core"
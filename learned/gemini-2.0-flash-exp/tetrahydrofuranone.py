"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is a tetrahydrofuran ring with a carbonyl group attached to any position of the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for THF ring pattern (C-C-C-C-O)
    thf_pattern = Chem.MolFromSmarts("C1CCOC1")
    if not mol.HasSubstructMatch(thf_pattern):
        return False, "No tetrahydrofuran ring found"

    # Check for carbonyl group attached to ring using updated SMARTS
    carbonyl_pattern = Chem.MolFromSmarts("[C,O]~[C](=O)")
    
    ring_atoms = mol.GetSubstructMatches(thf_pattern)[0]
    carbonyl_found = False
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if mol.GetSubstructMatch(Chem.MolFromSmarts(f"[{atom.GetSymbol()}]~[C](=[O])"), neighbor.GetIdx()):
                carbonyl_found = True
                break
        if carbonyl_found:
            break
    
    if not carbonyl_found:
         return False, "No carbonyl group attached to the THF ring"
    
    return True, "Contains a tetrahydrofuran ring with a carbonyl group attached to the ring"
"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:30736 tetrahydrofuranone
Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane with an oxo- substituent on the tetrahydrofuran ring.

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
    
    # Look for oxolane (tetrahydrofuran) ring with specific stereochemistry
    oxolane_pattern = Chem.MolFromSmarts("[O;R]1[C@@H]([C@@H](C)C)[C@@H](C)CC1")
    if not mol.HasSubstructMatch(oxolane_pattern):
        return False, "No oxolane (tetrahydrofuran) ring with correct stereochemistry found"
    
    # Look for oxo- substituent on the ring
    oxo_pattern = Chem.MolFromSmarts("[O;X1]=[C;R]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    
    # Check if any oxo- substituent is attached to the oxolane ring
    for oxo_match in oxo_matches:
        oxo_atom = mol.GetAtomWithIdx(oxo_match[1])
        for bond in oxo_atom.GetBonds():
            if bond.GetBeginAtomIdx() == oxo_match[1]:
                ring_atom_idx = bond.GetEndAtomIdx()
            else:
                ring_atom_idx = bond.GetBeginAtomIdx()
            
            ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
            if ring_atom.IsInRingSize(5) and ring_atom.HasProp("_ChiralityPossible"):
                return True, "Contains oxolane (tetrahydrofuran) ring with oxo- substituent"
    
    # Additional filter: molecular weight should be between 100 and 500 Da
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, "Molecular weight outside the expected range for tetrahydrofuranones"
    
    return False, "No oxo- substituent attached to the oxolane (tetrahydrofuran) ring"
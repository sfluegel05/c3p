"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone is defined as any isoflavonoid with a 3-aryl-1-benzopyran-4-one 
    (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoflavone core SMARTS pattern
    isoflavone_smarts = '[O]=C1C=CC2=CC=CC=C2O1'  # Chromen-4-one core
    isoflavone_core = Chem.MolFromSmarts(isoflavone_smarts)
    if isoflavone_core is None:
        return False, "Error in isoflavone core SMARTS pattern"

    # Check if the molecule contains the chromen-4-one core
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "Does not contain chromen-4-one core structure"

    # Define SMARTS pattern for aryl group attached at position 3
    # We will label atom indices in the core to find the attachment point
    isoflavone_core_labeled = Chem.MolFromSmarts('[O]=C1[C;H]=C[C;H]=C2O1[C;H]=C[C;H]=C2')  # Labeled chromen-4-one core
    matches = mol.GetSubstructMatches(isoflavone_core_labeled)
    if not matches:
        return False, "Chromen-4-one core not found with correct labeling"

    # For each match, check if there is an aryl group attached at position 3
    aryl_smarts = '[a]'

    for match in matches:
        # Get the atom at position 3 (index 1 in SMARTS pattern)
        pos3_atom_idx = match[1]
        pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
        # Get neighbors of position 3 atom
        neighbors = pos3_atom.GetNeighbors()
        aryl_found = False
        for neighbor in neighbors:
            if neighbor.GetIdx() not in match:
                # Check if neighbor is part of an aryl group
                neighbor_atom = neighbor.GetIdx()
                path = Chem.rdmolops.GetShortestPath(mol, pos3_atom_idx, neighbor_atom)
                submol = Chem.PathToSubmol(mol, path)
                if submol.HasSubstructMatch(Chem.MolFromSmarts(aryl_smarts)):
                    aryl_found = True
                    break
        if aryl_found:
            return True, "Contains isoflavone core with 3-aryl substitution"
    
    return False, "No 3-aryl substitution found on chromen-4-one core"
"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids with a specific cucurbitane framework.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    # Attempt to parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated SMARTS pattern for recognizing the cucurbitane backbone
    # This is a generic tetracyclic pattern; it may need further refinement based on known cucurbitane structures
    cucurbitane_pattern = Chem.MolFromSmarts('C[C@]1([C@@H]2C[C@@H]3[C@H]4C[C@@H](C)C5C=C(O)C[C@@]5(C)C4C=CC=C3C2C[C@]1(C)C(=O)CO')
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "No tetracyclic cucurbitane backbone structure found"

    # Checking for at least two keto groups (C=O)
    keto_pattern = Chem.MolFromSmarts('[CX3]=O')
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    if len(keto_matches) < 2:
        return False, f"Insufficient keto groups, found {len(keto_matches)}"

    # Cucurbitacins usually have a few hydroxyls (-OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"

    # Verify though less strict on molecular weight this time
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a typical cucurbitacin"

    return True, "Contains cucurbitane tetracyclic backbone with characteristic functional groups"
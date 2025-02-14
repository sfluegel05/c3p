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
    
    # Corrected SMARTS pattern for recognizing the cucurbitane backbone, which is C27H41
    cucurbitane_pattern = Chem.MolFromSmarts('C1[C@H]2CC[C@@H]3C[C@@]4(C)C=CC5(C)CCC(O)C=C5[C@@H](O)C4C=C3C2[C@@]1(C)C')
    if not mol.HasSubstructMatch(cucurbitane_pattern):
        return False, "No tetracyclic cucurbitane backbone structure found"

    # Checking for at least two keto groups (C=O) using a more precise pattern
    keto_pattern = Chem.MolFromSmarts('[CX3](=O)')
    keto_matches = mol.GetSubstructMatches(keto_pattern)
    if len(keto_matches) < 2:
        return False, f"Insufficient keto groups, found {len(keto_matches)}"

    # Cucurbitacins usually have a few hydroxyls (-OH)
    # RDKit does not directly support detecting hydroxyls from SMILES, need to check directly from atom reads.
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    if hydroxyl_count < 2:
        return False, f"Insufficient hydroxyl groups, found {hydroxyl_count}"

    # Check molecular weight - cucurbitacins typically are larger due to multiple rings and substitutions
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a typical cucurbitacin"

    return True, "Contains cucurbitane tetracyclic backbone with characteristic functional groups"
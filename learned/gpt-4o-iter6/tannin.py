"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds often forming complex glucosides with catechol and pyrogallol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyphenolic structure: multiple aromatic rings with hydroxyls
    # Phenolic group pattern (para-hydroxyl phenol as basic unit of catechol/pyrogallol)
    phenol_pattern = Chem.MolFromSmarts("c1(O)cccc1O")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) < 2:
        return False, "Not enough phenolic groups (need at least 2)"

    # Glucoside linkage pattern: O-glycosidic bond
    glycoside_pattern = Chem.MolFromSmarts("COC(=O)C")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside link detected"

    # Count hydroxyl groups to establish polyphenolic nature
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Too few oxygens (hydroxyl groups), found {o_count}"

    # Optional: Check for ester linkages (catechol and pyrogallol-based esters)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]") 
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester linkage found (optional criterion)"

    return True, "Polyphenolic structure with sufficient hydroxyl and glycoside linkages detected"
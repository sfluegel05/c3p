"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    
    A tannin is identified by the presence of specific polyphenolic structures such as catechol and pyrogallol groups.

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

    # Define catechol (1,2-dihydroxybenzene) pattern
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    
    # Define pyrogallol (1,2,3-trihydroxybenzene) pattern
    pyrogallol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    
    # Check for presence of catechol or pyrogallol groups
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    pyrogallol_matches = mol.GetSubstructMatches(pyrogallol_pattern)
    
    if not catechol_matches and not pyrogallol_matches:
        return False, "No catechol or pyrogallol groups found"
    
    # Check for glucoside complexity indicator
    # Glucose substructure can be approximated by a 6-membered ring with oxygens
    glucose_pattern = Chem.MolFromSmarts("O1[C@H]([C@@H](O)[C@H](O)[C@H](O)[C@H]1O)")
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    
    if not glucose_matches:
        return False, "No glucose or similar glucoside structures found"

    # Additional complexity checks could involve molecular weight or number of rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, "Not enough aromatic rings for tannin"

    return True, "Contains polyphenolic structures (catechol/pyrogallol) indicative of tannins with glucoside complexity"

# Add metadata and additional characterization if needed
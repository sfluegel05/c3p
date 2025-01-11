"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds with multiple hydroxyl groups on aromatic rings,
    often with glycosidic or ester linkages or galloyl groups.

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

    # Count aromatic rings
    aromatic_rings = sum(1 for ring in Chem.GetSymmSSSR(mol) if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if aromatic_rings < 1:
        return False, "No aromatic rings found"

    # Count phenolic hydroxyl groups (-OH attached to aromatic rings)
    phenolic_oh_pattern = Chem.MolFromSmarts("[OH][c]")
    phenolic_oh_matches = mol.GetSubstructMatches(phenolic_oh_pattern)
    if len(phenolic_oh_matches) < 2:
        return False, f"Found {len(phenolic_oh_matches)} phenolic hydroxyl groups, need at least 2"

    # Check for polyphenolic structure (multiple aromatic rings with hydroxyl groups)
    polyphenolic_pattern = Chem.MolFromSmarts("[c][OH].[c][OH]")
    polyphenolic_matches = mol.GetSubstructMatches(polyphenolic_pattern)
    if len(polyphenolic_matches) < 1:
        return False, "No polyphenolic structure found"

    # Check molecular weight - tannins are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for tannin"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for tannin"
    if o_count < 4:
        return False, "Too few oxygens for tannin"

    return True, "Contains multiple phenolic hydroxyl groups on aromatic rings, indicating a polyphenolic structure typical of tannins"
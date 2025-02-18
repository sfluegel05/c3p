"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group (OP(O)(=O)O or OP(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for polyprenol chain (at least 3 isoprene units)
    # More flexible pattern to account for different bond types and stereochemistry
    isoprene_pattern = Chem.MolFromSmarts("[C]=[C]-[C]-[C]-[C]=[C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 3"

    # Check that the phosphate is attached to the terminal carbon of the polyprenol chain
    # The phosphate should be connected to a carbon with a single hydroxyl group
    # More specific pattern to ensure proper attachment
    terminal_phosphate_pattern = Chem.MolFromSmarts("[C]([OX2]P(=O)([OX2])[OX2])~[C]=[C]")
    if not mol.HasSubstructMatch(terminal_phosphate_pattern):
        return False, "Phosphate not properly attached to terminal carbon of polyprenol chain"

    # Check for long carbon chain (at least 10 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"

    # Check molecular weight - polyprenol phosphates typically have MW > 300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for polyprenol phosphate"

    # Additional check for terminal double bond near phosphate
    terminal_double_bond_pattern = Chem.MolFromSmarts("[C]=[C]~[C]([OX2]P(=O)([OX2])[OX2])")
    if not mol.HasSubstructMatch(terminal_double_bond_pattern):
        return False, "No terminal double bond near phosphate group"

    return True, "Contains polyprenol chain with terminal phosphate group properly attached"
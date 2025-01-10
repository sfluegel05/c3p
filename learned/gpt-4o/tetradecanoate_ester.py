"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester
    based on its SMILES string.
    A tetradecanoate ester is obtained by condensation of
    myristic acid (tetradecanoic acid) with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the tetradecanoic acid moiety pattern using SMILES: A 14-carbon chain ending in a carboxyl
    tetradecanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    # Define ester linkage pattern, ensuring it connects to the carboxyl group
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[Cc,!#6;R1]") # Ensures attachment to a non-carbon single bonded atom and identifies ring structures

    # Check for tetradecanoic acid moiety
    if not mol.HasSubstructMatch(tetradecanoic_acid_pattern):
        return False, "No tetradecanoic acid moiety found"

    # Find all matches for tetradecanoic acid moiety
    myristic_matches = mol.GetSubstructMatches(tetradecanoic_acid_pattern)

    # Check for ester linkage that connects to a carboxyl group
    for match in myristic_matches:
        # Check each match to ensure ester linkage is correctly attached to this specific moiety
        myristic_atom_idx = match[-1]  # Get the index of the oxygen atom in the carboxyl
        for atom in mol.GetAtomWithIdx(myristic_atom_idx).GetNeighbors():
            neighbor_idx = atom.GetIdx()
            neighbor_smarts = f"[*:1][C:2](=O)O[*:3]"
            neighbor_pattern = Chem.MolFromSmarts(neighbor_smarts)
            if mol.GetSubstructMatch(neighbor_pattern):
                return True, "Contains tetradecanoic acid moiety with ester linkage"

    return False, "Ester linkage not matched with tetradecanoic acid moiety"
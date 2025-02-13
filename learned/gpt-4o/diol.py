"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol contains exactly two hydroxy groups, generally assumed part of its main chain.
    Complex cases with more hydroxys should have identifiable separate roles for two groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find hydroxyl groups (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Check primarily that there are two significant hydroxy groups as a diol
    if len(hydroxy_matches) == 2:
        return True, "Contains exactly two hydroxy groups, qualifying it as a diol"
    
    # Complex logic: check if adjacent diol groups or other structural cues
    elif len(hydroxy_matches) > 2:
        # Check if any two OH are adjacent or structurally relevant
        adjacent_oxygens = False
        # Create bond pattern queries
        o_bond_query = rdqueries.BondBetweenAtomsQuery()
        for i, j in hydroxy_matches:
            # For more than two OH groups, ensure at least two are connected or near-main chain
            for bond in mol.GetBonds():
                if (bond.GetBeginAtomIdx() in (i, j) and
                   bond.GetEndAtomIdx() in (i, j)):
                    adjacent_oxygens = True
                    break
            if adjacent_oxygens:
                break
        if adjacent_oxygens:
            return True, "Contains adjacent two hydroxy groups or structurally significant arrangement"
    
    return False, f"Found {len(hydroxy_matches)} hydroxy groups, which is inadequate for typical diol classification"
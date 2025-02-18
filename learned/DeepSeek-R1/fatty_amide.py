"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI fatty amide (monocarboxylic acid amide derived from a fatty acid)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid (aliphatic chain >=4 carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Find all amide groups (including primary amides)
    amide_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[NX3H2,H1,H0]')
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    # Check for monocarboxylic acid origin (only one carboxylic acid-derived amide)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2H1]')
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains carboxylic acid group"

    # Preferentially check primary amides first
    primary_amide_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[NH2]')
    primary_amide_matches = mol.GetSubstructMatches(primary_amide_pattern)
    
    for match in primary_amide_matches + amide_matches:
        carbonyl_idx = match[0]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Find the R group carbon attached to the carbonyl
        r_carbon = None
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in {match[1], match[2]}:
                r_carbon = neighbor
                break
        if not r_carbon:
            continue

        # Get longest aliphatic chain using RDKit's built-in function
        chain = rdMolDescriptors._CalcCarbons(mol, rootAtom=r_carbon.GetIdx(), onlyOnAromatic=False)
        chain_length = len(chain)
        
        # Check chain length and aliphatic nature
        if chain_length >= 4 and all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in chain):
            return True, f"Amide with {chain_length}-carbon aliphatic chain"

    return False, "No amide with sufficient aliphatic chain length"
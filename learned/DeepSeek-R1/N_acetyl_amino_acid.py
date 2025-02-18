"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:76295 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group (CH3-C=O) attached to the alpha-amino group
    of an amino acid (carboxylic acid group on adjacent carbon).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Verify exactly one carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(acid_matches) != 1:
        return False, f"Found {len(acid_matches)} carboxylic acid groups (needs exactly 1)"

    # Get the alpha carbon (carbon adjacent to carboxylic acid)
    acid_carbon = acid_matches[0][0]
    alpha_carbon = [n.GetIdx() for n in mol.GetAtomWithIdx(acid_carbon).GetNeighbors() 
                   if n.GetAtomicNum() == 6 and n.GetIdx() != acid_carbon]
    
    if not alpha_carbon:
        return False, "No alpha carbon adjacent to carboxylic acid"
    alpha_carbon = alpha_carbon[0]

    # Check for acetylated amine attached to alpha carbon
    # Pattern matches N-C(=O)-C (acetyl) where N is attached to alpha carbon
    acetyl_amide = Chem.MolFromSmarts("[N&H0;!$(N-C=O)]-&@[C&H0](-[C&H3])=O")
    for match in mol.GetSubstructMatches(acetyl_amide):
        nitrogen_idx = match[0]
        # Verify nitrogen is bonded to alpha carbon
        if mol.GetAtomWithIdx(nitrogen_idx).GetNeighbors()[0].GetIdx() == alpha_carbon:
            # Ensure acetyl group (exactly 1 methyl attached to carbonyl)
            carbonyl_carbon = match[1]
            methyl_count = sum(1 for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors()
                              if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalNumHs() >= 3)
            if methyl_count == 1:
                return True, "Acetylated amine on alpha carbon with carboxylic acid"

    return False, "No proper acetylated alpha-amino group found"
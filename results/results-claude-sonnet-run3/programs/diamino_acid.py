from rdkit import Chem
from rdkit.Chem import AllChem

def is_diamino_acid(smiles: str):
    """
    Determines if a molecule is a diamino acid (amino acid with two amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diamino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Count number of amino groups (including primary and secondary amines)
    amine_pattern = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    
    if len(amine_matches) < 2:
        return False, f"Only {len(amine_matches)} amino groups found - requires at least 2"
    
    # Check if amino groups are attached to carbon atoms
    valid_amino_groups = 0
    amino_positions = []
    
    for amine_idx in amine_matches:
        amine_atom = mol.GetAtomWithIdx(amine_idx[0])
        for neighbor in amine_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                valid_amino_groups += 1
                # Get position relative to carboxylic acid
                path = Chem.GetShortestPath(mol, neighbor.GetIdx(), 
                                          mol.GetSubstructMatch(carboxylic_pattern)[0])
                if path:
                    amino_positions.append(len(path))
                break
                
    if valid_amino_groups < 2:
        return False, "Less than 2 carbon-bound amino groups found"

    # Sort positions for consistent reporting
    amino_positions.sort()
    
    return True, f"Diamino acid with amino groups at positions {amino_positions} relative to carboxylic acid"
# Pr=0.9166666666666666
# Recall=1.0
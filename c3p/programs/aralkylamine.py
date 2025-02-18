"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) indicating if it's an aralkylamine and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find amine: primary, secondary or tertiary
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)]")
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No primary, secondary, or tertiary amine group found"

    # Find aromatic group
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAromaticAtoms()]
    if not aromatic_atoms:
        return False, "No aromatic group found"
    
    # Check connectivity between amine, alkyl chain, and aromatic
    aralkylamine_pattern = Chem.MolFromSmarts("[NX3][C][C][c]") # Simplified aryl-alkyl connection
    if not mol.HasSubstructMatch(aralkylamine_pattern):
        return False, "No alkyl group substituted by aromatic group found attached to amine."
   
    # Additional pattern to capture more complex structures, ensuring connectivity
    complex_aralkylamine_pattern = Chem.MolFromSmarts("[NX3][C,c][C,c][a]") 
    if not (mol.HasSubstructMatch(aralkylamine_pattern) or mol.HasSubstructMatch(complex_aralkylamine_pattern)):
        return False, "No alkyl group substituted by aromatic group found attached to amine."
    
    return True, "Molecule is an aralkylamine"

# The refined patterns aim to be inclusive enough to capture the diversity of aralkylamines in question without including unrelated compounds.
"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: ganglioside
"""

from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is composed of a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked on the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ceramide substructure
    # General ceramide pattern: amide bond connected to a sphingoid base (long-chain amino alcohol)
    ceramide_smarts = 'C(=O)N[C][C][C](O)'  # Amide bond connected to chain with hydroxyl group
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide substructure found"
    
    # Check for sialic acid residue
    # Simplified sialic acid pattern: a sugar ring with carboxylic acid group
    sialic_acid_smarts = 'C(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O'  # Sialic acid pattern
    sialic_acid_pattern = Chem.MolFromSmarts(sialic_acid_smarts)
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid residue found"
    
    # Check for glycosidic bonds (sugar chain)
    glycosidic_bond_smarts = '[C;R][O;!R][C;R]'  # Oxygen between two ring carbons (glycosidic bond)
    glycosidic_bond_pattern = Chem.MolFromSmarts(glycosidic_bond_smarts)
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bond_matches) < 2:
        return False, "Insufficient glycosidic linkages (sugar chain)"
    
    return True, "Contains ceramide backbone, sialic acid residue, and glycosidic linkages"
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
    # Ceramide general pattern: long-chain fatty acid linked via amide bond to a sphingoid base (long-chain amino alcohol)
    ceramide_smarts = 'C(=O)N[C@@H]([C;R0])[C;R0][C;R0](O)[C;R0]=[C;R0][C;R0]'  # Simplified ceramide pattern
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide substructure found"

    # Check for sialic acid substructure
    # Generalized sialic acid pattern: nine-carbon sugar with carboxylic acid and amino group
    sialic_acid_smarts = 'C(=O)[O;H1,-]C[C@@H](O)[C@H](O)[C@@H](O)[C@H](O[C@H]1[C@H](O)[C@H](O)[C@@H](NC=O)[C@@H]1O)CO'  # Simplified sialic acid pattern
    sialic_acid_pattern = Chem.MolFromSmarts(sialic_acid_smarts)
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid residue found"

    # Check for glycosidic linkage between ceramide and sugar chain
    # Look for glycosidic bond between sphingoid base and sugar
    glycosidic_linkage_smarts = 'O[C@H]([C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O)[C@@H]([C;R0])NC(=O)'  # Simplified glycosidic linkage
    glycosidic_linkage_pattern = Chem.MolFromSmarts(glycosidic_linkage_smarts)
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage between ceramide and sugar chain"

    return True, "Contains ceramide backbone, glycosidic-linked sugar chain, and sialic acid residue"
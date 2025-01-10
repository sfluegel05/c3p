"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetate ester SMARTS pattern: R-O-C(=O)C
    # Oxygens: [#8] with valence 2
    # Carbons: [#6] with appropriate bonding
    # The pattern represents an ester linkage where the carbonyl carbon is bonded to a methyl group
    acetate_ester_smarts = "[#8X2H0]-[#6](=O)-[CH3]"
    acetate_ester_pattern = Chem.MolFromSmarts(acetate_ester_smarts)

    if acetate_ester_pattern is None:
        return False, "Error in SMARTS pattern"

    # Search for the acetate ester pattern
    if mol.HasSubstructMatch(acetate_ester_pattern):
        return True, "Contains acetate ester group"
    else:
        return False, "No acetate ester group found"


__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'acetate ester',
        'definition': 'Any carboxylic ester where the carboxylic acid component is acetic acid.',
    },
    'examples': {
        'positive': [
            'CC(=O)OC[C@H]1O[C@@H](Oc2cc3C(=O)c4cc(O)ccc4C(=O)c3c(O)c2CO)[C@H](O)[C@@H](O)[C@@H]1O',
            'COC(=O)C[C@H]1[C@@]2(C)C[C@@]3(OC(C)=O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](O)C4=C1[C@H](OC(C)=O)[C@@]3(OC(C)=O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1',
            'C[C@H]1CC[C@@H]([C@@]2([C@H]([C@]3(C[C@]12C)C(=O)OCC3=C)OC(C)=O)[H])OC(=O)[C@@H](C)CC',
            'CCC(C)C(=O)O[C@H]1[C@H](O)C[C@@H]2[C@@](C)([C@@H]3C[C@H]4CCO[C@H]4O3)[C@H](C)C[C@H](OC(C)=O)[C@@]2(COC(C)=O)[C@@]11CO1',
            'C(C/C=C/COC(C)=O)C',
            # ... (additional examples)
        ],
        'negative': [
            # List of SMILES strings that are not acetate esters
        ],
    },
}
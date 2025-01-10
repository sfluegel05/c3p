"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside includes a ceramide backbone and oligosaccharide with one or more sialic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefine ceramide backbone pattern (sphingosine with amide linkage)
    ceramide_pattern = Chem.MolFromSmarts("NC(C=O)[C@@H](O)C(C)CCCCCCCCCCCCC")
    
    # Define sialic acid pattern (extended to capture variations)
    sialic_acid_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)C=O)C(O)[C@H](CO)N[C@@H](C)C(=O)")

    # Check for ceramide backbone
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for one or more sialic acid residues
    if len(mol.GetSubstructMatches(sialic_acid_pattern)) < 1:
        return False, "No sialic acid residues found"

    # Improved detection of glycosidic bonds (even in branched sugars)
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O[C@H](O)[C@@H]([C@H]1O)O)O")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkages found"

    return True, "Molecule is a ganglioside with ceramide backbone, sialic acid, and oligosaccharide structure"
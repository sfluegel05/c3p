"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define ceramide backbone pattern (long hydrophobic alkyl chain)
    ceramide_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\C=C\CCCCCCCCCCC")

    # Define sialic acid pattern (e.g., N-acetylneuraminic acid)
    sialic_acid_pattern = Chem.MolFromSmarts("O[C@@H](C=O)C(O)[C@H](CO)[C@H](O)[C@H](NC(C)=O)CO")

    # Check for ceramide backbone 
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for sialic acid(s)
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid residue found"

    # Check for glycosidic bonds (general sugar linkage)
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O[C@H](CO)[C@@H](C1O)O)")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkages found"

    return True, "Molecule is a ganglioside with ceramide backbone and sialic acid linked via glycosidic bonds"
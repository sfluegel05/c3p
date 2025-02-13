"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: CHEBI:36357 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for common biological molecule fragments
    amino_acid_pattern = Chem.MolFromSmarts("[N;H2,H1;!$(NC=O)]")
    nucleic_acid_pattern = Chem.MolFromSmarts("N1C=NC2=NC=NC=N12")
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCC(=O)O")
    carbohydrate_pattern = Chem.MolFromSmarts("OC[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO")
    peptide_pattern = Chem.MolFromSmarts("[N;H2,H1]C(=O)[C@H]([N;H2,H1])C(=O)")

    biological_fragments = [
        mol.HasSubstructMatch(amino_acid_pattern),
        mol.HasSubstructMatch(nucleic_acid_pattern),
        mol.HasSubstructMatch(fatty_acid_pattern),
        mol.HasSubstructMatch(carbohydrate_pattern),
        mol.HasSubstructMatch(peptide_pattern)
    ]

    # Check for at least 2 different types of biological fragments
    num_fragment_types = sum(1 for frag in biological_fragments if frag)
    if num_fragment_types < 2:
        return False, "Less than 2 types of biological fragments found"

    # Check for covalent linkage between fragments
    linker_pattern = Chem.MolFromSmarts("[N,O,S]-!@[N,O,S]")
    if not mol.HasSubstructMatch(linker_pattern):
        return False, "No covalent linkage found between biological fragments"

    # Additional checks or rules can be added here

    return True, "Contains at least 2 types of biological fragments covalently linked"
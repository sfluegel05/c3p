"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least 2 distinct biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded and refined SMARTS patterns for biological motifs or moieties
    peptide_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")   # Typical peptide bond
    nucleoside_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")  # Purine nucleobase structure
    coenzyme_pattern = Chem.MolFromSmarts("P(=O)(O)OCCN")      # CoA fragment
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")           # Thioester bond
    glycosidic_pattern = Chem.MolFromSmarts("[C@]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O1)") # Simple sugar ring
    sulfur_linkage_pattern = Chem.MolFromSmarts("S[CX3](=O)[O]") # Sulfate esters or similar linkages
    disulfide_pattern = Chem.MolFromSmarts("SS")               # typical disulfide bond
    phospholipid_pattern = Chem.MolFromSmarts("O=P(O)(O)OCC")  # Phospholipid head analog
        
    # Check for presence of distinct biological motifs
    patterns = [
        peptide_pattern, nucleoside_pattern, coenzyme_pattern,
        thioester_pattern, glycosidic_pattern, sulfur_linkage_pattern,
        disulfide_pattern, phospholipid_pattern
    ]
    
    distinct_biomotifs = sum(bool(mol.HasSubstructMatch(pattern)) for pattern in patterns)

    if distinct_biomotifs < 2:
        return False, "Less than two distinct biological motifs found"
    
    # Check molecular size to ensure complexity typically associated with bioconjugates
    mol_weight = Descriptors.MolWt(mol)
    if mol_weight < 500 and distinct_biomotifs < 3:
        return True, "Contains multiple distinct biological motifs, even with lower molecular weight"
    
    return True, "Contains multiple distinct biological motifs and linkage features, indicative of a bioconjugate"
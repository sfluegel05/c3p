"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least 2 biological molecules covalently linked together.

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

    # Refined and additional biological substructure patterns to enhance accuracy
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)][CX4][CX3](=O)[O,N,R]")
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H]([O,N])[C@H]1O")
    lipid_chain_pattern = Chem.MolFromSmarts("[CH3][CH2]{7,}")
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)[N,O]")

    # Covalent linkage patterns (extended to capture amide, thioester, ether, etc.)
    ester_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O][#6]")
    amide_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")
    thioester_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[S][#6]")
    ether_linkage_pattern = Chem.MolFromSmarts("[C][O][C]")
    disulfide_linkage_pattern = Chem.MolFromSmarts("[S][S]")

    # Biological patterns map
    patterns = {
        "amino acid": amino_acid_pattern,
        "sugar": sugar_pattern,
        "lipid chain": lipid_chain_pattern,
        "nucleotide": nucleotide_pattern,
    }

    # Covalent link patterns map
    linkage_patterns = {
        "ester linkage": ester_linkage_pattern,
        "amide linkage": amide_linkage_pattern,
        "thioester linkage": thioester_linkage_pattern,
        "ether linkage": ether_linkage_pattern,
        "disulfide linkage": disulfide_linkage_pattern,
    }

    # Detected patterns
    detected_bio_patterns = set()
    detected_cov_link_patterns = set()

    # Detect biological structures
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            detected_bio_patterns.add(name)

    # Detect covalent links
    for name, pattern in linkage_patterns.items():
        if mol.HasSubstructMatch(pattern):
            detected_cov_link_patterns.add(name)

    # Determine bioconjugate status
    if len(detected_bio_patterns) >= 2 and len(detected_cov_link_patterns) > 0:
        return True, f"Contains {len(detected_bio_patterns)} biologically relevant substructures and {len(detected_cov_link_patterns)} types of covalent linkages"
    else:
        return False, f"Contains {len(detected_bio_patterns)} biologically relevant substructures and {len(detected_cov_link_patterns)} types of covalent linkages, need at least 2 biological substructures with at least 1 covalent linkage"
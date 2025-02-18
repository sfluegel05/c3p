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
        return None, "Invalid SMILES string"

    try:
        # Biological pattern SMARTS (e.g., amino acids, sugars, nucleotides)
        amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)][CX4][CX3](=O)[O,N,R]")
        sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H]([O,N])[C@H]1O")
        nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
        
        # Additional Biological patterns
        cofactor_pattern = Chem.MolFromSmarts("CSCCC(=O)O")
        ester_linkage_pattern = Chem.MolFromSmarts("[CX3](=O)[O][CX4]")
        
        # Covalent linkage patterns (e.g., amide, disulfide, thioester)
        amide_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")
        disulfide_linkage_pattern = Chem.MolFromSmarts("[S][S]")
        thioester_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")

        # Biological structures
        bio_patterns = {
            "amino acid": amino_acid_pattern,
            "sugar": sugar_pattern,
            "nucleotide": nucleotide_pattern,
            "cofactor": cofactor_pattern
        }

        # Linkage structures
        linkage_patterns = {
            "amide linkage": amide_linkage_pattern,
            "disulfide linkage": disulfide_linkage_pattern,
            "thioester linkage": thioester_linkage_pattern,
            "ester linkage": ester_linkage_pattern
        }

        # Check for biological substructures
        detected_bio_patterns = set()
        for name, pattern in bio_patterns.items():
            if pattern and mol.HasSubstructMatch(pattern):
                detected_bio_patterns.add(name)

        # Check for linkage types
        detected_link_patterns = set()
        for name, pattern in linkage_patterns.items():
            if pattern and mol.HasSubstructMatch(pattern):
                detected_link_patterns.add(name)

        # Determine bioconjugate status based on detected patterns
        if len(detected_bio_patterns) >= 2 and len(detected_link_patterns) > 0:
            return True, f"Contains {len(detected_bio_patterns)} biological structures and {len(detected_link_patterns)} covalent link types"
        else:
            return False, f"Contains {len(detected_bio_patterns)} biological structures and {len(detected_link_patterns)} covalent link types. Requires at least 2 biological substructures with 1 linkage."

    except Exception as e:
        return None, f"Error in pattern matching: {e}"
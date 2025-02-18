"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    try:
        # Define more comprehensive biological patterns (including more biologically relevant structures)
        amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)][CX4][CX3](=O)[O,N,R]")
        peptide_linkage_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[NX3]")
        sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](O)[C@H]([O,N])[C@H]1O")
        nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
        
        cofactor_patterns = [
            Chem.MolFromSmarts("CSCCC(=O)O"),  # General cofactor motif
            Chem.MolFromSmarts("COP(O)(=O)O")   # Phosphate-containing motif e.g., CoA
        ]
        
        # Covalent linkage patterns (e.g., amide, disulfide, thioester)
        amide_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")
        disulfide_linkage_pattern = Chem.MolFromSmarts("[S][S]")
        thioester_linkage_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")

        # Patterns for biological components
        bio_patterns = {
            "amino acid": amino_acid_pattern,
            "sugar": sugar_pattern,
            "nucleotide": nucleotide_pattern,
            "cofactor": cofactor_patterns
        }

        # Patterns for linkage types
        linkage_patterns = {
            "amide linkage": amide_linkage_pattern,
            "disulfide linkage": disulfide_linkage_pattern,
            "thioester linkage": thioester_linkage_pattern,
            "peptide linkage": peptide_linkage_pattern
        }

        # Checking detected biological components
        detected_bio_patterns = set()
        for name, pattern in bio_patterns.items():
            if isinstance(pattern, list):
                for p in pattern:
                    if p and mol.HasSubstructMatch(p):
                        detected_bio_patterns.add(name)
                        break
            elif pattern and mol.HasSubstructMatch(pattern):
                detected_bio_patterns.add(name)

        # Checking detected linkage types
        detected_link_patterns = set()
        for name, pattern in linkage_patterns.items():
            if pattern and mol.HasSubstructMatch(pattern):
                detected_link_patterns.add(name)

        # Decision: must have at least 2 biological components and 1 linkage type
        if len(detected_bio_patterns) >= 2 and len(detected_link_patterns) > 0:
            return True, f"Contains {len(detected_bio_patterns)} biological structures and {len(detected_link_patterns)} covalent link types"
        else:
            return False, f"Contains {len(detected_bio_patterns)} biological structures and {len(detected_link_patterns)} covalent link types. Requires at least 2 biological substructures with 1 linkage."

    except Exception as e:
        return None, f"Error in pattern matching: {e}"
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

    # Define a variety of biological substructure patterns
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2][CX4][CX3](=O)[OX2H1]")
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H]1O")
    lipid_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    peptide_bond_pattern = Chem.MolFromSmarts("N[C@@H](C)C=O")
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2ncnc(N)c12")

    # Find matches for these patterns
    patterns = {
        "amino acid": amino_acid_pattern,
        "sugar": sugar_pattern,
        "lipid chain": lipid_chain_pattern,
        "peptide bond": peptide_bond_pattern,
        "nucleotide": nucleotide_pattern,
    }

    # Create a set to keep track of detected patterns
    detected_patterns = set()

    for name, pattern in patterns.items():
        matches = mol.GetSubstructMatches(pattern)
        if len(matches) > 0:
            detected_patterns.add(name)

    # Check for common chemical bonds indicative of covalent linking
    covalent_link_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0][#6]")  # ester linkage
    if mol.HasSubstructMatch(covalent_link_pattern):
        detected_patterns.add("covalent link")

    # Check if at least 2 distinct biological patterns and covalent linkage are found
    if len(detected_patterns) >= 2 and "covalent link" in detected_patterns:
        return True, f"Contains {len(detected_patterns)} biologically relevant substructures covalently linked"
    else:
        return False, f"Contains {len(detected_patterns)} biologically relevant substructures, need at least 2 with covalent linkage"

# Example usage:
result, reason = is_bioconjugate("S(SS(O)(=O)=O)C[C@H](N)C(O)=O")
print(result, reason)
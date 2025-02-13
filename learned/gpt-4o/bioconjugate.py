"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two distinct biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a bioconjugate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define substructure patterns
    patterns = {
        "amino_acid": Chem.MolFromSmarts("[NX3][CX4](C)C(=O)[OX2H1]"),  # Improved amino acid pattern
        "fatty_acid": Chem.MolFromSmarts("C(=O)O[CX4][CX4]C(C)..."),    # Multi-bond to capture longer chains
        "nucleotide": Chem.MolFromSmarts("n1cnc2c1ncnc2"),             # DNA/RNA base
        "sugar": Chem.MolFromSmarts("C(O)C(O)C(O)C(O)"),                # Generic sugar
        "peptide_bond": Chem.MolFromSmarts("C(=O)N[C]"),                # Broad peptide bond angle
        "thiol": Chem.MolFromSmarts("C[S]"),                            # Generic thiol group
        "phosphate": Chem.MolFromSmarts("P(=O)(O)O")                    # Phosphate backbone
    }

    # Track identified biological moieties
    identified_moieties = set()

    # Check for the presence of each biological moiety pattern
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            identified_moieties.add(name)

    # Check that they are covalently linked
    if len(identified_moieties) >= 2:
        for atom in mol.GetAtoms():
            if any(atom.HasSubstructMatch(pattern) for pattern in patterns.values()):
                if len(set(pattern for pattern in patterns if atom.HasSubstructMatch(pattern))) >= 2:
                    return True, f"Contains at least two different biological moieties covalently linked: {', '.join(identified_moieties)}"

    return False, "Does not contain at least two distinct biological components"

# Example use (replace SMILES with an actual string to test):
# result, reason = is_bioconjugate("your_smiles_here")
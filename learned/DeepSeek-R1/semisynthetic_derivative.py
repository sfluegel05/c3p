"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: semisynthetic derivative (CHEBI:134730)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on SMILES.
    Criteria: Contains natural product core scaffold AND synthetic modifications.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Define common natural product scaffolds using SMARTS
    NATURAL_SCAFFOLDS = [
        # Beta-lactam antibiotic core
        "[C]1([N])[C][C]=O",  # Simplified beta-lactam
        # Macrolide core (large lactone ring, at least 12 atoms)
        "[O]C(=O)C1CCCCC[C@H](O)CCCC",  # Approximate macrolide pattern
        # Steroid skeleton
        "C1C2CC(C3C4CCC5CCCCC5C4C3CC2)C1",
        # Alkaloid-like (e.g., indole)
        "c1ccc2c(c1)cnc3ccccc23",
        # Tetracycline-like
        "C1C(=O)C2C(C3C(C(C(C3(C(=O)O)O)O)O)O)C(C(=O)N(C)C)C1(O)O2",
    ]

    # Check for any natural product scaffold
    has_natural_core = False
    natural_reason = ""
    for idx, smarts in enumerate(NATURAL_SCAFFOLDS):
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_natural_core = True
            natural_reason = ["Beta-lactam", "Macrolide", "Steroid", "Alkaloid", "Tetracycline"][idx]
            break

    if not has_natural_core:
        return False, "No natural product core detected"

    # Check for synthetic modifications
    SYNTHETIC_FEATURES = [
        ("[Cl,Br,F,I]", "halogen atom"),
        ("[N+](=O)[O-]", "nitro group"),
        ("C#N", "cyano group"),
        ("[S](=O)(=O)", "sulfonyl group"),
        ("C(F)(F)F", "trifluoromethyl group"),
    ]

    synthetic_features = []
    for smarts, desc in SYNTHETIC_FEATURES:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            synthetic_features.append(desc)

    if not synthetic_features:
        return False, f"Natural {natural_reason} core but no synthetic modifications"

    return True, f"Contains {natural_reason} core with synthetic features: {', '.join(synthetic_features)}"
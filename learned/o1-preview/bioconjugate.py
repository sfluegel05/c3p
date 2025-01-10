"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
"""
from rdkit import Chem

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

    # Define SMARTS patterns for common biological molecules
    patterns = {
        'amino acid': Chem.MolFromSmarts('N[C@@H](C)C(=O)O'),  # Amino acid backbone
        'peptide bond': Chem.MolFromSmarts('N[C@@H](C)C(=O)N[C@@H](C)C(=O)'),  # Peptide linkage
        'nucleotide': Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # Purine base
        'sugar': Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'),  # Sugar ring
        'fatty acid': Chem.MolFromSmarts('C(=O)OCCCC'),  # Fatty acid chain
        'glutathione': Chem.MolFromSmarts('NCC(=O)N[C@H](CS)C(=O)NCC(=O)O'),  # Glutathione tripeptide
        'coenzyme A': Chem.MolFromSmarts('C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc23)C(O)C1O'),  # CoA structure
        'lipid': Chem.MolFromSmarts('C(=O)O[CX4][CX4]'),  # Ester linkage in lipids
        'polyketide': Chem.MolFromSmarts('C(=O)[CX3][CX2][CX3](=O)'),  # Polyketide chain
        'sterol': Chem.MolFromSmarts('C1CCC2C(C1)CCC3C(C2)CCCC3'),  # Steroid backbone
    }

    # Initialize a set to store detected biological molecules
    bio_molecules = set()

    # Search for patterns in the molecule
    for name, pattern in patterns.items():
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            bio_molecules.add(name)

    # Count the number of different biological molecules found
    num_bio_molecules = len(bio_molecules)

    # Check if at least two different biological molecules are present
    if num_bio_molecules >= 2:
        reason = f"Contains at least two biological molecules covalently linked: {', '.join(sorted(bio_molecules))}"
        return True, reason
    else:
        reason = "Does not contain at least two distinct biological molecules covalently linked"
        if num_bio_molecules == 1:
            reason += f": only found {next(iter(bio_molecules))}"
        return False, reason

# Example usage:
# smiles = "C(CCCC)[C@@H]([C@@H](\\C=C\\[C@H](C/C=C\\C/C=C\\CCCC(O)=O)O)SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O)O"
# result = is_bioconjugate(smiles)
# print(result)
"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
        'amino acid': Chem.MolFromSmarts('N[C@@H](C)C(=O)O'),
        'peptide': Chem.MolFromSmarts('N[C@@H](C)C(=O)N[C@@H](C)C(=O)'),
        'nucleotide': Chem.MolFromSmarts('n1cnc2c1ncnc2'),
        'sugar': Chem.MolFromSmarts('C1(C(C(C(O1)O)O)O)O'),
        'fatty acid': Chem.MolFromSmarts('C(=O)OCCCCCCCC'),
        'glutathione': Chem.MolFromSmarts('NCC(=O)N[C@@H](CS)C(=O)NCC(=O)O'),
        'coenzyme A': Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O'),
        'lipid': Chem.MolFromSmarts('C(=O)O[C@@H]('),
        'polyketide': Chem.MolFromSmarts('C(=O)C(C(=O))'),
    }

    # Initialize a list to store detected biological molecules
    bio_molecules = []

    # Search for patterns
    for name, pattern in patterns.items():
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            bio_molecules.append(name)

    # Remove duplicates
    bio_molecules = list(set(bio_molecules))

    # Count the number of different biological molecules found
    num_bio_molecules = len(bio_molecules)

    # Check if at least two different biological molecules are present
    if num_bio_molecules >= 2:
        reason = f"Contains at least two biological molecules covalently linked: {', '.join(bio_molecules)}"
        return True, reason
    else:
        return False, "Does not contain at least two distinct biological molecules covalently linked"

# Example usage:
# smiles = "C(CCCC)[C@@H]([C@@H](\\C=C\\[C@H](C/C=C\\C/C=C\\CCCC(O)=O)O)SC[C@H](NC(CC[C@H](N)C(=O)O)=O)C(=O)NCC(=O)O)O"
# result = is_bioconjugate(smiles)
# print(result)
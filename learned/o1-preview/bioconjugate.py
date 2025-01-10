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

    # RDKit requires hydrogens to be explicit for certain operations
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for common biological molecules
    patterns = {
        'amino acid': Chem.MolFromSmarts('[NX3,NX4+][CX4;H1,H2][CX3](=O)[O-]'),  # General amino acid pattern
        'nucleotide': Chem.MolFromSmarts('n1cnc2c1ncnc2'),  # Purine base
        'sugar': Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O'),  # Generic sugar ring
        'fatty acid': Chem.MolFromSmarts('C(=O)OCCCCCCCC'),  # Fatty acid chain
        'glutathione': Chem.MolFromSmarts('NCC(=O)N[C@H](CS)C(=O)NCC(=O)O'),  # Glutathione tripeptide
        'coenzyme A': Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc23)C(O)C1O'),  # CoA
        'lipid': Chem.MolFromSmarts('C(=O)O[CX4][CX4]'),  # Ester linkage in lipids
        'polyketide': Chem.MolFromSmarts('[#6][CX3](=O)[#6][CX3](=O)[#6]'),  # Polyketide chain
        'peptide': Chem.MolFromSmarts('N[C@@H](C(=O))C'),  # Peptide bond
    }

    # Fragment the molecule using BRICS fragmentation
    try:
        frags = Chem.BRICSDecompose(mol)
    except:
        return False, "Error during fragmentation"

    # Initialize a set to store detected biological molecules
    bio_molecules = set()

    # Convert fragment SMILES back to molecules
    frag_mols = [Chem.MolFromSmiles(frag) for frag in frags if Chem.MolFromSmiles(frag) is not None]

    # Search for patterns in fragments
    for frag in frag_mols:
        for name, pattern in patterns.items():
            if pattern is None:
                continue
            if frag.HasSubstructMatch(pattern):
                bio_molecules.add(name)

    # Count the number of different biological molecules found
    num_bio_molecules = len(bio_molecules)

    # Check if at least two different biological molecules are present
    if num_bio_molecules >= 2:
        reason = f"Contains at least two biological molecules covalently linked: {', '.join(sorted(bio_molecules))}"
        return True, reason
    else:
        return False, "Does not contain at least two distinct biological molecules covalently linked"

# Example usage:
# smiles = "C(CCCC)[C@@H]([C@@H](\\C=C\\[C@H](C/C=C\\C/C=C\\CCCC(O)=O)O)SC[C@H](NC(=O)CC[C@H](N)C(=O)O)C(=O)NCC(=O)O)O"
# result = is_bioconjugate(smiles)
# print(result)
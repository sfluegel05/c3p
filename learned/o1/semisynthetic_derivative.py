"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is an organic molecule derived from a natural product by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is organic (contains carbon atoms)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Molecule is not an organic compound"

    # Generate Murcko Scaffold (core structure of molecule)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold, isomericSmiles=True)

    # Define known natural product scaffolds (simplified for this example)
    # In practice, this would be a more comprehensive list
    natural_product_scaffolds = [
        'C1CCC(CC1)',  # Cyclohexane ring (terpenoids)
        'C1=CC=CC=C1',  # Benzene ring (common in alkaloids)
        'C1CC2CC3CC1CC(C2)C3',  # Steroid nucleus
        'C1CC[C@H]2[C@@H]3CCC(C3)CC[C@]12C',  # Steroid nucleus with stereochemistry
        'O=C1CC[C@H]2[C@@H]3CC[C@@H](O)[C@H](O)[C@@H](CC3)CC[C@]12C',  # Example triterpene scaffold
    ]
    # Check if scaffold matches any known natural product scaffold
    has_natural_product_core = False
    for np_smiles in natural_product_scaffolds:
        np_mol = Chem.MolFromSmiles(np_smiles)
        if scaffold.HasSubstructMatch(np_mol):
            has_natural_product_core = True
            break

    if not has_natural_product_core:
        return False, "Molecule does not contain a known natural product core"

    # Check for synthetic modifications
    synthetic_modifications = False

    # List of functional groups that might indicate synthetic modification
    # Focus on modifications uncommon in natural biosynthesis
    synthetic_patterns = {
        'Halogenation': '[F,Cl,Br,I]',
        'Azide': 'N=[N+]=[N-]',
        'Alkyne': 'C#C',
        'Nitro': '[$([NX3](=O)=O)]',
        'Trifluoromethyl': 'C(F)(F)F',
        'Methoxy': 'OC',
        'Ester': 'C(=O)OC',
        'Amide': 'C(=O)N',  # Amides can be natural, but check for unusual patterns
        'Propargyl': 'C#CC',
        'Protecting group': 'C(=O)O[C@H1]',  # Simplified pattern
    }

    for mod_name, smarts in synthetic_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            synthetic_modifications = True
            break  # One modification is enough to consider

    if not synthetic_modifications:
        return False, "No synthetic modifications detected"

    # Ensure that the molecule is not a known natural product (i.e., unmodified)
    # We can compare the scaffold to the entire molecule
    if scaffold_smiles == Chem.MolToSmiles(mol, isomericSmiles=True):
        return False, "Molecule appears to be a natural product without modifications"

    return True, "Molecule contains a natural product core with evidence of synthetic modifications"
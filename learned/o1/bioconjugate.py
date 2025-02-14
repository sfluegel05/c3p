"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is a molecular entity consisting of at least 2 different biological molecules covalently linked together.

    Biological molecule classes considered:
    - Protein (includes amino acids and peptides)
    - Carbohydrate
    - Nucleotide
    - Lipid
    - Cofactor (e.g., Coenzyme A)

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

    types_found = set()

    # Define SMARTS patterns for biological molecule classes
    biomolecule_patterns = {
        'protein': [
            # Pattern to detect peptide bonds (amide bonds between amino acids)
            Chem.MolFromSmarts("N[C@@H](C)C(=O)N[C@@H](C)C(=O)"),  # Simplified peptide chain
            Chem.MolFromSmarts("N[C@@H](C)C(=O)O"),  # Amino acid residue
            Chem.MolFromSmarts("N[C@@H](C)C(=O)N")  # General peptide bond
        ],
        'carbohydrate': [
            # Pattern to detect pyranose and furanose rings with multiple hydroxyl groups
            Chem.MolFromSmarts("[$([C@@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)]"),  # Pyranose
            Chem.MolFromSmarts("[$([C@@H]1O[C@H](O)[C@@H](O)[C@H]1O)]"),  # Furanose
        ],
        'nucleotide': [
            # Pattern to detect nucleobases attached to sugar and phosphate
            Chem.MolFromSmarts("n1c[nH]c2c1ncnc2"),  # Purine base
            Chem.MolFromSmarts("c1[nH]cnc1"),  # Pyrimidine base
            Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@@H](CO)O1"),  # Ribose sugar
            Chem.MolFromSmarts("P(=O)(O)[O-]")  # Phosphate group
        ],
        'lipid': [
            # Pattern to detect long aliphatic chains (fatty acids)
            Chem.MolFromSmarts("C(=O)@[O;H1]@[C;X4][C;X4][C;X4][C;X4][C;X4]"),  # Fatty acid tail
            Chem.MolFromSmarts("C(=O)SCCN")  # Acyl-CoA linkage
        ],
        'cofactor': [
            # Pattern to detect Coenzyme A fragments
            Chem.MolFromSmarts("C1=NC2=C(N1)N=CN2C3C(C(C(O3)COP(=O)(O)O)O)O"),  # Adenine nucleotide
            Chem.MolFromSmarts("OP(=O)(O)OCCNC(=O)CCNC(=O)")  # Pantetheine moiety in CoA
        ]
    }

    # Check for presence of each biological molecule class
    for biomolecule_class, patterns in biomolecule_patterns.items():
        for pattern in patterns:
            if mol.HasSubstructMatch(pattern):
                types_found.add(biomolecule_class)
                break  # No need to check other patterns for this class

    # Determine if the molecule is a bioconjugate
    if len(types_found) >= 2:
        return True, f"Contains {', '.join(types_found)} covalently linked"
    else:
        return False, "Does not contain at least 2 different covalently linked biological units"
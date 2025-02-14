"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    The B vitamins include vitamin B1 (thiamine), B2 (riboflavin), B3 (niacin),
    B5 (pantothenic acid), B6 (pyridoxine, pyridoxal, pyridoxamine),
    B7 (biotin), B9 (folic acid), and B12 (cobalamin).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a B vitamin, False otherwise
        str: Reason for classification
    """

    # Convert input SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of known B vitamin SMILES and their names
    b_vitamin_smiles = [
        # Vitamin B1 (Thiamine)
        ('CC1=C(C)C=CN=C1CCO', 'Vitamin B1 (Thiamine)'),
        # Vitamin B2 (Riboflavin)
        ('CN(C)C1=NC2=C(N1)C(=O)N(C)C3=C2C=C(C=C3)O', 'Vitamin B2 (Riboflavin)'),
        # Vitamin B3 (Niacin)
        ('c1ccncc1C(=O)O', 'Vitamin B3 (Niacin - Nicotinic acid)'),
        ('c1ccncc1C(=O)N', 'Vitamin B3 (Niacin - Nicotinamide)'),
        # Vitamin B5 (Pantothenic acid)
        ('CC(C)(CO)C(=O)NCCC(=O)O', 'Vitamin B5 (Pantothenic acid)'),
        # Vitamin B6 (Pyridoxine, Pyridoxal, Pyridoxamine)
        ('O=Cc1ncc(CO)c(C)c1O', 'Vitamin B6 (Pyridoxal)'),
        ('NCc1ncc(CO)c(C)c1O', 'Vitamin B6 (Pyridoxamine)'),
        ('COc1ncc(CO)c(C)c1O', 'Vitamin B6 (Pyridoxine)'),
        # Vitamin B7 (Biotin)
        ('O=C1NC(=O)N2C[C@@H](SC1)[C@]2([H])CCCCCC(=O)O', 'Vitamin B7 (Biotin)'),
        # Vitamin B9 (Folic acid and derivatives)
        ('Nc1nc2ncc(CNc3ccc(cc3)C(=O)O)nc2c(=O)[nH]1', 'Vitamin B9 (Folic acid)'),
        # Vitamin B12 (Cobalamin)
        # Cobalamin is complex; check for cobalt atom coordinated in a corrin ring
        # Here we use a simplified pattern for cobalamin
        ('[Cobalt]', 'Vitamin B12 (Cobalamin)'),
    ]

    # Convert known B vitamin SMILES to molecules
    b_vitamin_mols = []
    for smi, name in b_vitamin_smiles:
        mol_vit = Chem.MolFromSmiles(smi)
        if mol_vit:
            b_vitamin_mols.append((mol_vit, name))

    # Check if input molecule matches any known B vitamin
    for vit_mol, vit_name in b_vitamin_mols:
        if vit_name == 'Vitamin B12 (Cobalamin)':
            # Check for cobalt atom
            contains_cobalt = any(atom.GetAtomicNum() == 27 for atom in mol.GetAtoms())
            if contains_cobalt:
                return True, f"Molecule matches {vit_name} (Contains cobalt atom characteristic of cobalamin)"
            else:
                continue
        else:
            # Perform substructure matching
            if mol.HasSubstructMatch(vit_mol):
                return True, f"Molecule matches {vit_name}"

    return False, "Molecule does not match any known B vitamin"
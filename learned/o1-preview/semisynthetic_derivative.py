"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:90773 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdFMCS import FindMCS

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

    # Generate Murcko scaffold of the molecule
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return False, "Unable to generate scaffold"

    # Find the Maximum Common Substructure (MCS) between the molecule and its scaffold
    mcs = FindMCS([mol, scaffold], matchValences=True, ringMatchesRingOnly=True)
    if mcs.numAtoms == 0:
        return False, "No common substructure found"

    # Calculate the proportion of the molecule that is modified
    mol_atoms = mol.GetNumHeavyAtoms()
    scaffold_atoms = scaffold.GetNumHeavyAtoms()
    mcs_atoms = mcs.numAtoms
    modified_atoms = mol_atoms - mcs_atoms

    # If a significant portion of the molecule is modified, consider synthetic modification
    if modified_atoms >= 5 and (modified_atoms / mol_atoms) >= 0.1:
        # Define functional groups commonly introduced synthetically
        synthetic_groups = [
            Chem.MolFromSmarts('[#6][F,Cl,Br,I]'),    # Halogenated carbons
            Chem.MolFromSmarts('c[F,Cl,Br,I]'),       # Halogenated aromatics
            Chem.MolFromSmarts('[#6][C#N]'),          # Nitriles
            Chem.MolFromSmarts('N(=O)[O-]'),          # Nitro group
            Chem.MolFromSmarts('S(=O)(=O)[O-]'),      # Sulfonic acids
            Chem.MolFromSmarts('C(=O)O[C,N]'),        # Esters and carbamates
            Chem.MolFromSmarts('C(=O)N[C,N]'),        # Amides
            Chem.MolFromSmarts('C(=O)O[C@H]'),        # Esters at chiral centers
            Chem.MolFromSmarts('[C;!R]=[C;!R]'),      # Unconjugated double bonds
            Chem.MolFromSmarts('O[C@@H1][#6]'),       # Ether linkages
            Chem.MolFromSmarts('C#C'),                # Alkynes
            Chem.MolFromSmarts('C=O'),                # Carbonyl groups
            Chem.MolFromSmarts('C[O,N,S]C'),          # Alkylation (ethers, amines)
            Chem.MolFromSmarts('C=C-C=O'),            # Enones
            Chem.MolFromSmarts('[#16R]'),             # Sulfur-containing rings
        ]

        # Check for synthetic modifications
        has_synthetic_modification = False
        modifications = []

        for group in synthetic_groups:
            if mol.HasSubstructMatch(group):
                has_synthetic_modification = True
                smarts = Chem.MolToSmarts(group)
                modifications.append(smarts)

        if has_synthetic_modification:
            reason = f"Contains synthetic modifications: {', '.join(modifications)}"
            return True, reason
        else:
            return False, "No synthetic functional groups detected"

    else:
        return False, "No significant synthetic modifications detected"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:90773',
        'name': 'semisynthetic derivative',
        'definition': 'Any organic molecular entity derived from a natural product by partial chemical synthesis.',
        'parents': []
    }
}
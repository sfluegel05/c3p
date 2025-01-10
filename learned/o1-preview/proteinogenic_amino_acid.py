"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    A proteinogenic amino acid is defined as:
    'Any of the 23 alpha-amino acids that are precursors to proteins, and are incorporated into proteins during translation.
    The group includes the 20 amino acids encoded by the nuclear genes of eukaryotes together with selenocysteine,
    pyrrolysine, and N-formylmethionine. Apart from glycine, which is non-chiral, all have L configuration.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Preserve isotopes and deuterium atoms
    # Remove explicit hydrogens to simplify matching
    mol = Chem.RemoveHs(mol)

    # Assign stereochemistry
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Generate canonical SMILES for input molecule
    input_canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    # Define canonical SMILES for all 23 proteinogenic amino acids
    amino_acids = {
        'Glycine': 'NCC(=O)O',
        'Alanine': 'N[C@@H](C)C(=O)O',
        'Valine': 'N[C@@H](C(C)C)C(=O)O',
        'Leucine': 'N[C@@H](CC(C)C)C(=O)O',
        'Isoleucine': 'N[C@@H](C(C)CC)C(=O)O',
        'Proline': 'N1CCCC1C(=O)O',
        'Phenylalanine': 'N[C@@H](CC1=CC=CC=C1)C(=O)O',
        'Tyrosine': 'N[C@@H](CC1=CC=C(O)C=C1)C(=O)O',
        'Tryptophan': 'N[C@@H](CC1=CNC2=CC=CC=C12)C(=O)O',
        'Serine': 'N[C@@H](CO)C(=O)O',
        'Threonine': 'N[C@@H](C(O)C)C(=O)O',
        'Cysteine': 'N[C@@H](CS)C(=O)O',
        'Methionine': 'N[C@@H](CCSC)C(=O)O',
        'Asparagine': 'N[C@@H](CC(=O)N)C(=O)O',
        'Glutamine': 'N[C@@H](CCC(=O)N)C(=O)O',
        'Lysine': 'N[C@@H](CCCCN)C(=O)O',
        'Arginine': 'N[C@@H](CCCNC(N)=N)C(=O)O',
        'Histidine': 'N[C@@H](CC1=CN=CN1)C(=O)O',
        'Aspartic Acid': 'N[C@@H](CC(=O)O)C(=O)O',
        'Glutamic Acid': 'N[C@@H](CCC(=O)O)C(=O)O',
        'Selenocysteine': 'N[C@@H](C[SeH])C(=O)O',
        'Pyrrolysine': 'N[C@@H](CCCCNC(=O)C1=CC=CN1)C(=O)O',
        'N-Formylmethionine': 'O=CN[C@@H](CCSC)C(=O)O',
    }

    # Convert the amino acid SMILES to RDKit molecules and prepare for matching
    amino_acid_mols = {}
    for name, aa_smiles in amino_acids.items():
        aa_mol = Chem.MolFromSmiles(aa_smiles)
        if aa_mol is not None:
            Chem.AssignAtomChiralTagsFromStructure(aa_mol)
            Chem.AssignStereochemistry(aa_mol, force=True, cleanIt=True)
            amino_acid_mols[name] = aa_mol

    # Compare input molecule to each amino acid
    for name, aa_mol in amino_acid_mols.items():
        # Create molecule copies to prevent modification
        mol_copy = Chem.Mol(mol)
        aa_mol_copy = Chem.Mol(aa_mol)
        # Use RDKit's isomorphic matching with chirality consideration
        if mol_copy.HasSubstructMatch(aa_mol_copy, useChirality=True):
            return True, f"Molecule matches {name}"

    return False, "Molecule is not a proteinogenic amino acid with L-configuration"
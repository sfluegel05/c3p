"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

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

    # Remove isotopic labels (e.g., deuterium, tritium) by setting isotope to zero
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)

    # Assign stereochemistry to ensure consistent chirality representation
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Generate canonical SMILES with stereochemistry
    input_canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    # Define canonical SMILES for all 23 proteinogenic amino acids (L-configuration)
    amino_acids = {
        'Glycine': 'NCC(=O)O',
        'Alanine': 'N[C@H](C)C(=O)O',
        'Valine': 'N[C@H](C(C)C)C(=O)O',
        'Leucine': 'N[C@H](CC(C)C)C(=O)O',
        'Isoleucine': 'N[C@H](C(C)CC)C(=O)O',
        'Proline': 'N1CCC[C@H]1C(=O)O',
        'Phenylalanine': 'N[C@H](CC1=CC=CC=C1)C(=O)O',
        'Tyrosine': 'N[C@H](CC1=CC=C(O)C=C1)C(=O)O',
        'Tryptophan': 'N[C@H](CC1=CNC2=CC=CC=C12)C(=O)O',
        'Serine': 'N[C@H](CO)C(=O)O',
        'Threonine': 'N[C@H](C(O)C)C(=O)O',
        'Cysteine': 'N[C@H](CS)C(=O)O',
        'Methionine': 'N[C@H](CCSC)C(=O)O',
        'Asparagine': 'N[C@H](CC(=O)N)C(=O)O',
        'Glutamine': 'N[C@H](CCC(=O)N)C(=O)O',
        'Lysine': 'N[C@H](CCCCN)C(=O)O',
        'Arginine': 'N[C@H](CCCNC(=N)N)C(=O)O',
        'Histidine': 'N[C@H](CC1=CN=C-N1)C(=O)O',
        'Aspartic Acid': 'N[C@H](CC(=O)O)C(=O)O',
        'Glutamic Acid': 'N[C@H](CCC(=O)O)C(=O)O',
        'Selenocysteine': 'N[C@H](C[SeH])C(=O)O',
        'Pyrrolysine': 'N[C@H](CCCCNC(=O)[C@H]1CCCN1)C(=O)O',
        'N-Formylmethionine': 'O=CN[C@H](CCSC)C(=O)O',
    }

    # Normalize and process each amino acid for comparison
    for name, aa_smiles in amino_acids.items():
        aa_mol = Chem.MolFromSmiles(aa_smiles)
        if aa_mol is not None:
            # Remove isotopes
            for atom in aa_mol.GetAtoms():
                atom.SetIsotope(0)
            # Assign stereochemistry
            Chem.AssignAtomChiralTagsFromStructure(aa_mol)
            Chem.AssignStereochemistry(aa_mol, force=True, cleanIt=True)
            # Generate canonical SMILES
            aa_canonical_smiles = Chem.MolToSmiles(aa_mol, isomericSmiles=True, canonical=True)
            # Compare canonical SMILES strings
            if input_canonical_smiles == aa_canonical_smiles:
                return True, f"Molecule matches {name}"

    return False, "Molecule is not a proteinogenic amino acid with L-configuration"
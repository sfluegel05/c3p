"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

# List of SMILES strings for the 23 proteinogenic amino acids
proteinogenic_amino_acids = [
    "N[C@@H](CCCNC(N)=N)C(O)=O",  # L-arginine
    "OC(=O)[C@@H](N)CC1=C(C(=C(O)C(=C1[2H])[2H])[2H])[2H]",  # L-tyrosine-d4
    "OC(=O)[C@@H](N)CC(C([2H])([2H])[2H])C",  # L-leucine-d3
    "N[C@@H](CO)C(O)=O",  # L-serine
    "N[C@@H](CS)C(O)=O",  # L-cysteine
    "N[C@@H](CCC(O)=O)C(O)=O",  # L-glutamic acid
    "N[C@@H](Cc1c[nH]cn1)C(O)=O",  # L-histidine
    "CSCC[C@H](N)C(O)=O",  # L-methionine
    "OC(=O)[C@]1(NC(C(C1([2H])[2H])([2H])[2H])([2H])[2H])[2H]",  # L-proline-d7
    "C(C(N([2H])[2H])([2H])[2H])(=O)O[2H]",  # glycine-d5
    "N[C@@H](CC(O)=O)C(O)=O",  # L-aspartic acid
    "N[C@@H](CC(N)=O)C(O)=O",  # L-asparagine
    "OC([C@H]([C@H](CC)C)N)=O",  # L-isoleucine
    "CC(C)[C@H](N)C(O)=O",  # L-valine
    "O[13C](=O)[13CH2][15NH2]",  # glycine-13C2,15N
    "O(C(=O)[C@@](N([2H])[2H])(C(C(O[2H])=O)([2H])[2H])[2H])[2H]",  # L-aspartic acid-d7
    "C[C@@H](O)[C@H](N)C(O)=O",  # L-threonine
    "NCCCC[C@H](N)C(O)=O",  # L-lysine
    "CC(C)C[C@H](N)C(O)=O",  # L-leucine
    "OC(=O)[C@@H]1CCCN1",  # L-proline
    "S(CC[C@H](N)C(O)=O)C([2H])([2H])[2H]",  # L-methionine-d3
    "C[C@H](N)C(O)=O",  # L-alanine
    "OC(=O)[C@@H](N)CC1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H]",  # L-phenylalanine-d5
    "OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]",  # L-arginine-d7
    "OC(=O)[C@@](N)(C(C([2H])([2H])[2H])(C([2H])([2H])[2H])[2H])[2H]",  # L-valine-d8
    "N[C@@H](CCC(N)=O)C(O)=O",  # L-glutamine
    "C(=O)([C@@H](N)CCCCNC([C@H]1[C@@H](CC=N1)C)=O)O"  # L-pyrrolysine
]

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino and carboxyl groups
    has_amino_group = any(atom.GetSymbol() == 'N' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 for atom in mol.GetAtoms())
    has_carboxyl_group = any(atom.GetSymbol() == 'C' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 3 and sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == 'O') == 2 for atom in mol.GetAtoms())
    if not (has_amino_group and has_carboxyl_group):
        return False, "Missing amino and/or carboxyl group"

    # Check for chiral center with L/S configuration
    chiral_centers = Chem.FindMolChiralUnspecifiedUnbondedAtoms(mol)
    if not chiral_centers:
        return False, "No chiral center found"
    chiral_center = chiral_centers[0]
    if mol.GetAtomWithIdx(chiral_center).GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        return False, "Chiral center does not have L/S configuration"

    # Check if the molecule matches a known proteinogenic amino acid
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    if canonical_smiles in proteinogenic_amino_acids:
        return True, f"Matches proteinogenic amino acid: {canonical_smiles}"

    return False, "Molecule does not match any known proteinogenic amino acid"
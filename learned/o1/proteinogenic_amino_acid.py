"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:59870 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    A proteinogenic amino acid is one of the 23 alpha-amino acids that are incorporated into proteins during translation.
    Apart from glycine, which is non-chiral, all have L configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to ensure proper valence
    mol = Chem.AddHs(mol)

    # Assign stereochemistry
    Chem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define amino acid backbone SMARTS patterns
    # General pattern for alpha-amino acids (including isotopes and labels)
    backbone_smarts = '[N;X3H2,+0][C@?H](*)C(=O)[O;H1,-1]'
    backbone = Chem.MolFromSmarts(backbone_smarts)

    # Define glycine pattern (non-chiral)
    glycine_smarts = '[N;X3H2,+0][CH2]C(=O)[O;H1,-1]'
    glycine = Chem.MolFromSmarts(glycine_smarts)

    # Check for glycine first
    if mol.HasSubstructMatch(glycine):
        # Ensure the molecule is glycine and not part of a larger structure
        if mol.GetNumAtoms() > 7:
            return False, "Molecule is larger than glycine"
        return True, "Molecule is glycine, a proteinogenic amino acid"

    # Check for alpha-amino acid backbone
    matches = mol.GetSubstructMatches(backbone)
    if matches:
        # Ensure the molecule is a single amino acid (no peptide bonds)
        amide_bonds = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                if (atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6) or \
                   (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 7):
                    if atom1.GetHybridization() == Chem.rdchem.HybridizationType.SP2 or \
                       atom2.GetHybridization() == Chem.rdchem.HybridizationType.SP2:
                        amide_bonds += 1
        if amide_bonds > 1:
            return False, "Molecule contains peptide bonds, not a single amino acid"

        # Identify the alpha carbon atom
        for match in matches:
            alpha_carbon_idx = match[1]
            alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

            # Check chirality of alpha carbon
            chiral_tag = alpha_carbon.GetChiralTag()
            if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW or chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                # Get the stereochemistry (R or S)
                stereochemistry = Chem.rdMolDescriptors.CalcChiralCenters(mol, force=True, includeUnassigned=True)
                for center in stereochemistry:
                    idx, config = center
                    if idx == alpha_carbon_idx:
                        side_chain_atom_symbols = [nbr.GetSymbol() for nbr in alpha_carbon.GetNeighbors() if nbr.GetAtomicNum() > 1]
                        # For cysteine and selenocysteine, L-configuration corresponds to R
                        if 'S' in side_chain_atom_symbols or 'Se' in side_chain_atom_symbols:
                            if config == 'R':
                                return True, "Molecule has L-configuration at alpha carbon, a proteinogenic amino acid"
                            else:
                                return False, "Cysteine and selenocysteine must have R configuration at alpha carbon"
                        else:
                            if config == 'S':
                                return True, "Molecule has L-configuration at alpha carbon, a proteinogenic amino acid"
                            else:
                                return False, "Amino acids must have L-configuration at alpha carbon"
                return False, "Could not determine chirality"
            else:
                return False, "Alpha carbon is not properly chiral"
    else:
        return False, "Molecule does not have alpha-amino acid backbone"

    return False, "Molecule does not match proteinogenic amino acid"
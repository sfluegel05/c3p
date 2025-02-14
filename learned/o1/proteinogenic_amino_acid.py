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

    # Define alpha-amino acid backbone SMARTS pattern with variable side chain (R)
    backbone_smarts = '[NX3][C@H](*)C(=O)[OH]'
    backbone = Chem.MolFromSmarts(backbone_smarts)

    # Check if molecule matches the backbone (excluding glycine)
    matches = mol.GetSubstructMatches(backbone)
    if matches:
        # For chiral amino acids, check that the alpha carbon has L-configuration
        match = matches[0]
        alpha_carbon_idx = match[1]
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)

        # Check chirality of alpha carbon
        chiral_tag = alpha_carbon.GetChiralTag()
        if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW or chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            # Assign stereochemistry
            Chem.AssignAtomChiralTagsFromStructure(mol)
            stereochemistry = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
            for center in stereochemistry:
                idx, config = center
                if idx == alpha_carbon_idx:
                    # For most amino acids, L-configuration corresponds to S
                    # For cysteine, L-configuration corresponds to R due to sulfur priority
                    side_chain_atoms = [nbr for nbr in alpha_carbon.GetNeighbors() if nbr.GetIdx() not in match]
                    if len(side_chain_atoms) != 1:
                        return False, "Could not identify side chain"
                    side_chain_atom = side_chain_atoms[0]
                    if side_chain_atom.GetSymbol() == 'S':
                        # Cysteine case
                        if config == 'R':
                            return True, "Molecule is L-cysteine, a proteinogenic amino acid"
                        else:
                            return False, "Cysteine must have R configuration at alpha carbon"
                    else:
                        if config == 'S':
                            return True, "Molecule has L-configuration at alpha carbon, a proteinogenic amino acid"
                        else:
                            return False, "Amino acids must have L-configuration at alpha carbon"
            return False, "Could not determine chirality"
        else:
            return False, "Alpha carbon is not properly chiral"
    else:
        # Check for glycine, which is achiral and has two hydrogens on the alpha carbon
        glycine_smarts = '[NX3][CH2]C(=O)[OH]'
        glycine = Chem.MolFromSmarts(glycine_smarts)
        glycine_matches = mol.GetSubstructMatches(glycine)
        if glycine_matches:
            return True, "Molecule is glycine, a proteinogenic amino acid"
        else:
            return False, "Molecule does not have alpha-amino acid backbone"

    return False, "Molecule does not match proteinogenic amino acid"
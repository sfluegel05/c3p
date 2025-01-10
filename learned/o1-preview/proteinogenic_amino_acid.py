"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Normalize isotopes and convert deuterium to hydrogen
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            atom.SetIsotope(0)
        if atom.GetAtomicNum() == 1 and atom.GetMass() > 1.0079:
            # Convert deuterium ([2H]) to hydrogen
            atom.SetAtomicNum(1)
            atom.SetIsotope(0)
            atom.SetFormalCharge(0)

    # Remove explicit hydrogens to simplify matching
    mol = Chem.RemoveHs(mol)

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define the alpha-amino acid backbone SMARTS pattern
    # Matches alpha carbon connected to:
    # - Amino group [N;H2]
    # - Carboxyl group C(=O)O
    # - Side chain [*]
    alpha_amino_acid_pattern = Chem.MolFromSmarts('[N;H2][C@H](*)C(=O)O')

    # Check for presence of the alpha-amino acid backbone with L-configuration
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if matches:
        # Verify that alpha carbon has the correct stereochemistry
        for match in matches:
            alpha_c_idx = match[1]
            alpha_c = mol.GetAtomWithIdx(alpha_c_idx)
            if alpha_c.HasProp('_CIPCode'):
                cip_code = alpha_c.GetProp('_CIPCode')
                # For most amino acids, L-configuration corresponds to 'S' configuration
                # For cysteine and selenocysteine, L corresponds to 'R' due to sulfur/selenium priority
                neighbor_atomic_nums = [nbr.GetAtomicNum() for nbr in alpha_c.GetNeighbors()]
                if 16 in neighbor_atomic_nums or 34 in neighbor_atomic_nums:
                    # Contains sulfur (S, atomic number 16) or selenium (Se, atomic number 34)
                    if cip_code == 'R':
                        return True, "Molecule is L-cysteine or L-selenocysteine"
                else:
                    if cip_code == 'S':
                        return True, "Molecule is a proteinogenic amino acid with L-configuration"
    else:
        # Check for glycine (achiral)
        glycine_pattern = Chem.MolFromSmarts('NCC(=O)O')
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "Molecule is glycine (achiral proteinogenic amino acid)"
        # Check for N-formylmethionine and pyrrolysine explicitly
        n_formylmethionine = Chem.MolFromSmiles('O=CNCC[C@H](N)C(=O)O')
        pyrrolysine = Chem.MolFromSmiles('O=C(O)[C@@H](N)CCCCNC(=O)C1=CC=CN1')
        if mol.HasSubstructMatch(n_formylmethionine):
            return True, "Molecule is N-formylmethionine"
        if mol.HasSubstructMatch(pyrrolysine):
            return True, "Molecule is pyrrolysine"

    return False, "Molecule is not a proteinogenic amino acid with L-configuration"
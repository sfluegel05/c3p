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
    Apart from glycine, which is non-chiral, all have L configuration. This function checks for the standard backbone and side chains.

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

    # Add hydrogens (including isotopes) to ensure proper valence
    mol = Chem.AddHs(mol, addCoords=True)

    # Assign stereochemistry
    Chem.AssignAtomChiralTagsFromStructure(mol)
    AllChem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define the alpha-amino acid backbone pattern
    # Include primary amine [N;!H0] and secondary amine in case of proline [N;R]
    backbone_smarts = '[N;!H0]C[C;!H0](C(=O)[O])'
    backbone = Chem.MolFromSmarts(backbone_smarts)

    if not mol.HasSubstructMatch(backbone):
        return False, "Molecule does not have alpha-amino acid backbone"

    # Check for glycine (non-chiral)
    glycine_smarts = 'NCC(=O)O'  # Glycine backbone
    glycine = Chem.MolFromSmarts(glycine_smarts)
    if mol.HasSubstructMatch(glycine):
        return True, "Molecule is glycine, a proteinogenic amino acid"

    # Get chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found; amino acid (except glycine) must be chiral"

    # Get the alpha carbon (the one bonded to both amine and carboxylate groups)
    alpha_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon atom
            # Check if carbon is connected to nitrogen and carboxylate
            has_nitrogen = False
            has_carboxylate = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7:
                    has_nitrogen = True
                elif neighbor.GetAtomicNum() == 6:
                    # Check if this neighbor is part of a carboxylate group
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetAtomicNum() == 8:
                            has_carboxylate = True
            if has_nitrogen and has_carboxylate:
                alpha_carbons.append(atom)

    if not alpha_carbons:
        return False, "Could not find alpha carbon"

    # Assume the first alpha carbon found is the one we're interested in
    alpha_carbon = alpha_carbons[0]
    alpha_carb_idx = alpha_carbon.GetIdx()

    # Get the chiral tag of the alpha carbon
    chiral_tag = alpha_carbon.GetChiralTag()
    if chiral_tag == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        return False, "Alpha carbon is not chiral"

    # Determine if the configuration is correct
    # For L-amino acids (except cysteine and selenocysteine), the configuration is S
    # For cysteine and selenocysteine, the configuration is R

    # Determine the side chain (R group) attached to the alpha carbon
    side_chain_atoms = []
    for neighbor in alpha_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() != 1 and neighbor.GetAtomicNum() not in [7, 6]:
            side_chain_atoms.append(neighbor)

    if not side_chain_atoms:
        # Proline's side chain loops back to the nitrogen
        side_chain_atoms = [neighbor for neighbor in alpha_carbon.GetNeighbors() if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_carbon.GetIdx()]

    if not side_chain_atoms:
        return False, "No side chain found"

    # Create a Mol object of the side chain
    side_chain_idxs = Chem.GetMolFrags(Chem.PathToSubmol(mol, [alpha_carbon.GetIdx()] + [atom.GetIdx() for atom in side_chain_atoms]), asMols=False)[0]
    side_chain = Chem.PathToSubmol(mol, side_chain_idxs)

    # Define side chain SMARTS patterns for the 23 amino acids
    side_chain_smarts_list = {
        'Alanine': '[CH3]',
        'Arginine': 'CCCNC(N)=N',
        'Asparagine': 'CC(N)=O',
        'Aspartic acid': 'CC(=O)O',
        'Cysteine': 'CS[*]',  # Include isotopes
        'Glutamine': 'CCC(N)=O',
        'Glutamic acid': 'CCC(=O)O',
        'Glycine': '',  # Handled separately
        'Histidine': 'Cc1c[nH]cn1',
        'Isoleucine': 'C(C)CC',
        'Leucine': 'CC(C)C',
        'Lysine': 'CCCCN',
        'Methionine': 'CCSC',
        'Phenylalanine': 'CCc1ccccc1',
        'Proline': 'C1CCCN1',  # Proline's side chain forms a ring
        'Serine': 'CO',
        'Threonine': 'C(O)C',
        'Tryptophan': 'CCc1c[nH]c2ccccc12',
        'Tyrosine': 'CCc1ccc(O)cc1',
        'Valine': 'C(C)C',
        'Selenocysteine': 'C[Se]',
        'Pyrrolysine': 'CCCCNC(=O)C1CCCN1',  # Simplified side chain
        'N-formylmethionine': 'CCSC'  # Same side chain as methionine
    }

    # Check if side chain matches any of the standard amino acids
    side_chain_matched = False
    for name, pattern in side_chain_smarts_list.items():
        if pattern == '':
            continue  # Glycine handled before
        side_chain_mol = Chem.MolFromSmarts(pattern)
        if side_chain_mol is None:
            continue
        if side_chain.HasSubstructMatch(side_chain_mol):
            amino_acid_name = name
            side_chain_matched = True
            break

    if not side_chain_matched:
        return False, "Side chain does not match any of the 23 proteinogenic amino acids"

    # Determine the stereochemistry at the alpha carbon
    stereo_config = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    alpha_carbon_stereo = None
    for idx, config in stereo_config:
        if idx == alpha_carbon.GetIdx():
            alpha_carbon_stereo = config
            break

    if alpha_carbon_stereo is None:
        return False, "Could not determine chirality at alpha carbon"

    # For cysteine and selenocysteine, the L-configuration corresponds to R
    if amino_acid_name in ['Cysteine', 'Selenocysteine']:
        if alpha_carbon_stereo == 'R':
            return True, f"Molecule is {amino_acid_name} with L-configuration (R) at alpha carbon"
        else:
            return False, f"{amino_acid_name} must have R configuration at alpha carbon"
    else:
        # For all other amino acids, L-configuration corresponds to S
        if alpha_carbon_stereo == 'S':
            return True, f"Molecule is {amino_acid_name} with L-configuration (S) at alpha carbon"
        else:
            return False, f"{amino_acid_name} must have S configuration at alpha carbon"
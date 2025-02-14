"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:33709 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has an alpha-carbon with L-configuration (usually 'S' configuration, 'R' for cysteine/selenocysteine), attached to an amino group and a carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules with peptide bonds (amide linkages)
    peptide_bond = Chem.MolFromSmarts('C(=O)N')
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Molecule contains peptide bonds, likely a peptide not a single amino acid"

    # Add hydrogens to the molecule to ensure accurate stereochemistry
    mol = Chem.AddHs(mol)

    # Assign stereochemistry to the molecule
    AllChem.AssignAtomChiralTagsFromStructure(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Define a SMARTS pattern for the alpha-amino acid backbone
    # This pattern searches for a carbon (alpha-carbon) bonded to:
    # - an amino group [NH2 or NH]
    # - a carboxyl group [C(=O)OH]
    # - a side chain (any other atom)
    # The alpha-carbon should be chiral except for glycine
    alpha_amino_acid_pattern = Chem.MolFromSmarts('[N;H2,H1][C@@?H](*)C(=O)[O;H1,H0-]')
    if alpha_amino_acid_pattern is None:
        return False, "Error creating SMARTS pattern"

    # Search for matches of the alpha-amino acid backbone in the molecule
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if not matches:
        return False, "No alpha-amino acid backbone found"

    # Check that there is only one alpha-amino acid backbone
    if len(matches) > 1:
        return False, "Multiple alpha-amino acid backbones found"

    # Iterate over the matches
    for match in matches:
        # Get the index of the alpha-carbon atom (second atom in the SMARTS pattern)
        alpha_c_atom_idx = match[1]
        alpha_c_atom = mol.GetAtomWithIdx(alpha_c_atom_idx)

        # Check if the alpha-carbon is chiral (for all amino acids except glycine)
        if alpha_c_atom.HasProp('_CIPCode'):
            cip_code = alpha_c_atom.GetProp('_CIPCode')

            # Determine if the molecule is L-alpha-amino acid based on CIP code
            # Cysteine and selenocysteine are special cases (R configuration)
            side_chain_atoms = [nbr.GetAtomicNum() for nbr in alpha_c_atom.GetNeighbors() if nbr.GetIdx() not in [match[0], match[2]]]
            has_sulfur = 16 in side_chain_atoms  # 16 is sulfur's atomic number
            has_selenium = 34 in side_chain_atoms  # 34 is selenium's atomic number

            if cip_code == 'S':
                if has_sulfur or has_selenium:
                    return False, "Molecule is D-cysteine or D-selenocysteine (L-form has R configuration)"
                else:
                    return True, "Molecule is an L-alpha-amino acid with (S) configuration at alpha-carbon"
            elif cip_code == 'R':
                if has_sulfur or has_selenium:
                    return True, "Molecule is L-cysteine or L-selenocysteine with (R) configuration at alpha-carbon"
                else:
                    return False, "Molecule is D-alpha-amino acid"
            else:
                return False, "Unknown stereochemistry at alpha-carbon"
        else:
            # If alpha-carbon is not chiral, check if it's glycine (no side chain)
            neighbors = [nbr for nbr in alpha_c_atom.GetNeighbors() if nbr.GetIdx() not in [match[0], match[2]] and nbr.GetAtomicNum() > 1]
            if len(neighbors) == 0:
                # No side chain, likely glycine
                return True, "Molecule is glycine, an achiral alpha-amino acid"
            else:
                return False, "Alpha-carbon is not chiral or missing stereochemistry information"

    return False, "No appropriate alpha-amino acid backbone with correct stereochemistry found"
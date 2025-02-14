"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has an alpha-carbon with L-configuration, attached to an amino group and a carboxyl group.
    
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

    # Add hydrogens to the molecule to ensure accurate stereochemistry
    mol = Chem.AddHs(mol)

    # Assign stereochemistry to the molecule
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define a SMARTS pattern for the alpha-amino acid backbone with a chiral center
    alpha_amino_acid_pattern = Chem.MolFromSmarts('[N;!H0][C@H](*)C(=O)[O;H1,H0-]')
    if alpha_amino_acid_pattern is None:
        return False, "Error creating SMARTS pattern"

    # Search for matches of the alpha-amino acid backbone in the molecule
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if not matches:
        return False, "No alpha-amino acid backbone with chiral center found"

    # Iterate over all matches to check for L-configuration
    for match in matches:
        # Get the index of the alpha-carbon atom (second atom in the SMARTS pattern)
        alpha_c_atom_idx = match[1]
        alpha_c_atom = mol.GetAtomWithIdx(alpha_c_atom_idx)

        # Check if the alpha-carbon has CIP stereochemistry assigned
        if alpha_c_atom.HasProp('_CIPCode'):
            cip_code = alpha_c_atom.GetProp('_CIPCode')

            # Get neighboring atoms to the alpha-carbon (excluding amino nitrogen and carboxyl carbon)
            neighbors = [nbr for nbr in alpha_c_atom.GetNeighbors() if nbr.GetIdx() not in [match[0], match[2]] and nbr.GetAtomicNum() != 1]
            
            # There should be exactly one side chain atom (R group)
            if len(neighbors) != 1:
                continue  # Unable to determine side chain, skip to next match
            side_chain_atom = neighbors[0]

            # Check if the side chain contains sulfur (for cysteine and selenocysteine)
            has_sulfur = side_chain_atom.GetAtomicNum() == 16 or any(nb.GetAtomicNum() == 16 for nb in side_chain_atom.GetNeighbors())

            # Determine if the molecule is L-alpha-amino acid based on CIP code and side chain
            if cip_code == 'S':
                if has_sulfur:
                    return False, "Molecule is D-cysteine (L-cysteine has R configuration)"
                else:
                    return True, "Molecule is an L-alpha-amino acid with (S) configuration at alpha-carbon"
            elif cip_code == 'R':
                if has_sulfur:
                    return True, "Molecule is L-cysteine or L-selenocysteine with (R) configuration at alpha-carbon"
                else:
                    return False, "Molecule is D-alpha-amino acid"
        else:
            return False, "Alpha-carbon is not chiral or missing stereochemistry information"

    return False, "No appropriate alpha-amino acid backbone with correct stereochemistry found"
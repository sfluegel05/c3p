"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule containing four amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplify the molecule by removing side chain variations
    try:
        # Attempt to find the peptide backbone
        backbone_pattern = Chem.MolFromSmarts("[N;!R][C;!R][C;!R](=O)")
        backbone_matches = mol.GetSubstructMatches(backbone_pattern)

        # Count the number of amino acid residues
        residues = []
        for match in backbone_matches:
            residue_atoms = set(match)
            # Extend the residue to include side chains
            atom = mol.GetAtomWithIdx(match[1])  # Alpha carbon
            residue = Chem.GetMolFragment(mol, atom.GetIdx(), asMols=False, sanitizeFrags=False)
            residues.append(residue)

        num_residues = len(backbone_matches)

        if num_residues != 4:
            return False, f"Found {num_residues} amino acid residues, expected 4 for tetrapeptide"

        # Ensure residues are connected in sequence
        # Check for continuous peptide chain without interruptions
        if not Chem.MolToSmiles(mol).count("C(=O)N") >= 3:
            return False, "Peptide bonds do not form a continuous chain"

        # Check for N-terminal and C-terminal groups (optional)
        # Allows for variations like acetylation or amidation
        return True, "Molecule is a tetrapeptide with 4 amino acid residues connected via peptide bonds"

    except Exception as e:
        return False, f"Error during analysis: {str(e)}"
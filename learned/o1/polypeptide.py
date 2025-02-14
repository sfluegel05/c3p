"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide - A peptide containing ten or more amino acid residues.
"""

from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide containing ten or more amino acid residues)
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for peptide backbone
    # Pattern for the backbone nitrogen atom in peptide bond
    peptide_nitrogen_smarts = "[NX3;H1;!$(NC=[#6,#7,#8]);$(N-C-C(=O))]"
    peptide_nitrogen = Chem.MolFromSmarts(peptide_nitrogen_smarts)

    # Pattern for alpha carbon connected to nitrogen and carbonyl carbon
    alpha_carbon_smarts = "[CX4H1;$(C-[#6,#7,#8]);$(C-N-C(=O))]"
    alpha_carbon = Chem.MolFromSmarts(alpha_carbon_smarts)

    # Find all nitrogen atoms that are part of the peptide backbone
    peptide_nitrogens = mol.GetSubstructMatches(peptide_nitrogen)
    peptide_nitrogen_indices = [match[0] for match in peptide_nitrogens]

    # Find all alpha carbon atoms
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon)
    alpha_carbon_indices = [match[0] for match in alpha_carbons]

    # Identify amino acid residues by finding alpha carbons connected to peptide nitrogens
    residue_count = 0
    for idx in alpha_carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = atom.GetNeighbors()
        has_peptide_nitrogen = False
        has_carbonyl_carbon = False
        for neighbor in neighbors:
            if neighbor.GetIdx() in peptide_nitrogen_indices:
                has_peptide_nitrogen = True
            elif neighbor.GetAtomicNum() == 6:
                # Check if this carbon is a carbonyl carbon
                for bond in mol.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        begin_atom = bond.GetBeginAtom()
                        end_atom = bond.GetEndAtom()
                        if (begin_atom.GetIdx() == neighbor.GetIdx() and end_atom.GetAtomicNum() == 8) or \
                           (end_atom.GetIdx() == neighbor.GetIdx() and begin_atom.GetAtomicNum() == 8):
                            has_carbonyl_carbon = True
                            break
        if has_peptide_nitrogen and has_carbonyl_carbon:
            residue_count +=1

    if residue_count >= 10:
        return True, f"Contains {residue_count} amino acid residues connected via peptide bonds"
    else:
        return False, f"Contains {residue_count} amino acid residues, less than 10 required for polypeptide"
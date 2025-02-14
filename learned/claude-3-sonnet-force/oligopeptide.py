"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:36322 oligopeptide
"A peptide containing a relatively small number of amino acids."
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Define patterns and lists for amino acids and modifications
AA_SMARTS = ['[N;!$(NC=O)]', '[$(N(-C(=O))-C)]']  # Amino group, amino acid backbone
AA_NAMES = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His',
            'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
            'Tyr', 'Val']
MODS = ['C(=O)N', 'C(=O)O', 'OC(=O)', 'NC(=O)']  # Common terminal modifications

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide backbone and amino group
    has_peptide_backbone = any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in AA_SMARTS)
    if not has_peptide_backbone:
        return False, "No peptide backbone or amino group found"

    # Count amino acid residues
    aa_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N' and
                   sum(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom.GetNeighbors()) == 0)

    # Check for common terminal modifications
    has_mods = any(mol.HasSubstructMatch(Chem.MolFromSmarts(s)) for s in MODS)

    # Check for amino acid side chains
    aa_sidechains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and not mol.GetAtomWithIdx(atom.GetIdx()).GetIsAromatic():
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, atom.GetIdx(), 4)
            atom_dict = {atom_idx: mol.GetAtomWithIdx(atom_idx).GetSymbol() for atom_idx in env}
            sidechain = ''.join(sorted(atom_dict.values()))
            if sidechain not in aa_sidechains:
                aa_sidechains.append(sidechain)

    # Classify as oligopeptide based on criteria
    if 2 <= aa_count <= 20 and len(aa_sidechains) >= 2 and has_mods:
        return True, "Contains 2-20 amino acid residues with peptide bonds and common modifications"
    elif 2 <= aa_count <= 20 and len(aa_sidechains) >= 2:
        return True, "Contains 2-20 amino acid residues with peptide bonds"
    else:
        return False, "Does not meet criteria for an oligopeptide"
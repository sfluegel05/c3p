"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Standard proteinogenic amino acids in their neutral form
    proteinogenic = {
        'ALA': 'N[CH](C)C(=O)O',  # Alanine
        'ARG': 'N[CH](CCCNC(=N)N)C(=O)O',  # Arginine
        'ASN': 'N[CH](CC(=O)N)C(=O)O',  # Asparagine
        'ASP': 'N[CH](CC(=O)O)C(=O)O',  # Aspartic acid
        'CYS': 'N[CH](CS)C(=O)O',  # Cysteine
        'GLN': 'N[CH](CCC(=O)N)C(=O)O',  # Glutamine
        'GLU': 'N[CH](CCC(=O)O)C(=O)O',  # Glutamic acid
        'GLY': 'NCC(=O)O',  # Glycine
        'HIS': 'N[CH](CC1=CN=CN1)C(=O)O',  # Histidine
        'ILE': 'N[CH]([CH](CC)C)C(=O)O',  # Isoleucine
        'LEU': 'N[CH](CC(C)C)C(=O)O',  # Leucine
        'LYS': 'N[CH](CCCCN)C(=O)O',  # Lysine
        'MET': 'N[CH](CCSC)C(=O)O',  # Methionine
        'PHE': 'N[CH](Cc1ccccc1)C(=O)O',  # Phenylalanine
        'PRO': 'N1[CH](CCC1)C(=O)O',  # Proline
        'SER': 'N[CH](CO)C(=O)O',  # Serine
        'THR': 'N[CH]([CH](O)C)C(=O)O',  # Threonine
        'TRP': 'N[CH](CC1=CNC2=C1C=CC=C2)C(=O)O',  # Tryptophan
        'TYR': 'N[CH](Cc1ccc(O)cc1)C(=O)O',  # Tyrosine
        'VAL': 'N[CH](C(C)C)C(=O)O'  # Valine
    }

    # Check for basic amino acid structure (more flexible patterns)
    patterns = [
        '[NX3,NX4;!$(NC=O)][CX4][CX3](=[OX1])[OX2H,OX1-,OX2]',  # Standard pattern
        '[NX3,NX4;!$(NC=O)][CX4][CX3](=[OX1])[OX2R]',  # Ester form
        '[NX3,NX4]([H])[CX4][CX3](=[OX1])[OX2H,OX1-]'  # Explicit H
    ]
    
    has_amino_acid_core = False
    for pattern in patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_amino_acid_core = True
            break
            
    if not has_amino_acid_core:
        return False, "Missing basic amino acid structure"

    # Convert input to canonical SMILES
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Check if it's a standard amino acid (allowing for different protonation states)
    for aa_smiles in proteinogenic.values():
        prot_mol = Chem.MolFromSmiles(aa_smiles)
        if prot_mol is None:
            continue
        # Compare without stereochemistry and hydrogen counts
        prot_canonical = Chem.MolToSmiles(prot_mol, isomericSmiles=False)
        mol_canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
        if prot_canonical == mol_canonical:
            return False, "This is a standard proteinogenic amino acid"
    
    # Additional checks
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 4:
        return False, "Molecule too small to be an amino acid"
    
    # Count carbons and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 2:
        return False, "Insufficient carbon atoms for amino acid structure"
    if n_count < 1:
        return False, "Must contain at least one nitrogen"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 75:
        return False, f"Molecular weight {mol_wt:.1f} too small for amino acid"
    if mol_wt > 2000:
        return False, f"Molecular weight {mol_wt:.1f} too large for typical amino acid"
    
    # Look for common modifications that indicate non-proteinogenic status
    modifications = [
        ('Phosphorylation', '[PX4](=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]'),
        ('Methylation', '[NX3,NX4;!$(NC=O)](C)(C)'),
        ('Halogenation', '[F,Cl,Br,I]'),
        ('Hydroxylation', '[OX2H]'),
        ('Unusual cyclization', '[R3,R4,R5,R6,R7,R8]')
    ]
    
    found_modifications = []
    for mod_name, mod_smarts in modifications:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(mod_smarts)):
            found_modifications.append(mod_name)
    
    reason = "Non-proteinogenic amino acid identified"
    if found_modifications:
        reason += f" with modifications: {', '.join(found_modifications)}"
    
    return True, reason
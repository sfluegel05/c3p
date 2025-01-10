"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: non-proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Standard proteinogenic amino acids SMILES (simplified forms)
    proteinogenic = {
        'G': 'NCC(=O)O',  # Glycine
        'A': 'CC(N)C(=O)O',  # Alanine
        'R': 'NC(=N)NCCCC(N)C(=O)O',  # Arginine
        'N': 'NC(=O)CC(N)C(=O)O',  # Asparagine
        'D': 'OC(=O)CC(N)C(=O)O',  # Aspartic acid
        'C': 'SCC(N)C(=O)O',  # Cysteine
        'Q': 'NC(=O)CCC(N)C(=O)O',  # Glutamine
        'E': 'OC(=O)CCC(N)C(=O)O',  # Glutamic acid
        'H': 'NC(Cc1[nH]cnc1)C(=O)O',  # Histidine
        'I': 'CC(C)C(C)C(N)C(=O)O',  # Isoleucine
        'L': 'CC(C)CC(N)C(=O)O',  # Leucine
        'K': 'NCCCC(N)C(=O)O',  # Lysine
        'M': 'CSCC(N)C(=O)O',  # Methionine
        'F': 'NC(Cc1ccccc1)C(=O)O',  # Phenylalanine
        'P': 'OC(=O)C1CCCN1',  # Proline
        'S': 'OCC(N)C(=O)O',  # Serine
        'T': 'CC(O)C(N)C(=O)O',  # Threonine
        'W': 'NC(Cc1c[nH]c2ccccc12)C(=O)O',  # Tryptophan
        'Y': 'NC(Cc1ccc(O)cc1)C(=O)O',  # Tyrosine
        'V': 'CC(C)C(N)C(=O)O'  # Valine
    }
    
    # Check for basic amino acid structure
    amino_acid_pattern = Chem.MolFromSmarts('[NX3,NX4;!$(NC=O)][CX4][CX3](=[OX1])[OX2H,OX1-]')
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Missing basic amino acid structure (NH2-CH-COOH)"
    
    # Convert input to canonical SMILES for comparison
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    
    # Check if it's one of the standard amino acids
    for aa_smiles in proteinogenic.values():
        prot_mol = Chem.MolFromSmiles(aa_smiles)
        if prot_mol is None:
            continue
        # Convert to canonical form without stereochemistry for base structure comparison
        prot_canonical = Chem.MolToSmiles(prot_mol, isomericSmiles=False)
        mol_canonical = Chem.MolToSmiles(mol, isomericSmiles=False)
        if prot_canonical == mol_canonical:
            return False, "This is a standard proteinogenic amino acid"
    
    # Additional checks for common non-proteinogenic features
    
    # Count number of atoms (exclude very small molecules)
    if mol.GetNumAtoms() < 4:
        return False, "Molecule too small to be an amino acid"
    
    # Verify presence of carbon backbone
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 2:
        return False, "Insufficient carbon atoms for amino acid structure"
    
    # Check for reasonable molecular weight (most amino acids are between 75-500 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 75 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for amino acids"

    return True, "Non-proteinogenic amino acid structure identified"
"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    Criteria include multiple amide bonds in peptide backbone structure and molecular features.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all amide bonds (peptide bonds)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H0,H1]")
    amide_count = len(mol.GetSubstructMatches(amide_pattern))
    if amide_count < 5:
        return False, f"Insufficient amide bonds ({amide_count} < 5)"
    
    # Check for alpha-amino acid patterns (CH-R-C(=O)-N)
    # This helps distinguish peptide backbones from other amides
    alpha_aa_pattern = Chem.MolFromSmarts("[CH1X4]([!H0])[CX3](=O)[NX3]")
    if not mol.HasSubstructMatch(alpha_aa_pattern):
        return False, "No alpha-amino acid residues detected"
    
    # Check molecular weight range (500-3000 Da typical for peptide antibiotics)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 3000:
        return False, f"Molecular weight {mol_wt:.1f} Da outside typical range (500-3000 Da)"
    
    # Check for presence of non-proteinogenic amino acids or modified groups
    # Common in antibiotics: D-amino acids, N-methylation, cyclic structures
    # Look for at least one N-methyl amide or non-alpha amino acid as complexity indicator
    modified_amide = Chem.MolFromSmarts("[NX3H0]([CH3])[CX3](=O)")  # N-methyl amide
    if mol.HasSubstructMatch(modified_amide):
        return True, "Contains modified peptide features (e.g. N-methylation)"
    
    # Check for cyclic structure through bridged amides
    sssr = Chem.GetSSSR(mol)
    for ring in sssr:
        ring_atoms = list(ring)
        amide_in_ring = any(mol.GetBondBetweenAtoms(i,j).GetBondType() == Chem.BondType.AMIDE 
                           for i in ring_atoms for j in ring_atoms if i < j)
        if amide_in_ring and len(ring_atoms) <= 10:
            return True, "Contains cyclic peptide structure"
    
    # Fallback check for extended peptide chain
    peptide_chain = Chem.MolFromSmarts("[NX3H0,H1][CH1X4][CX3](=O)")
    chain_matches = len(mol.GetSubstructMatches(peptide_chain))
    if chain_matches >= 3:
        return True, "Contains extended peptide backbone structure"
    
    return False, "Does not meet structural criteria for peptide antibiotics"
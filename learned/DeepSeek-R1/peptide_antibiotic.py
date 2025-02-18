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
    Criteria include presence of multiple amide bonds and molecular weight.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amide bond pattern (peptide bonds)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H0,H1]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 4:
        return False, f"Only {len(amide_matches)} amide bonds found (needs at least 4)"
    
    # Check molecular weight (peptide antibiotics are typically >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    # Optional: Check for presence of amino and carboxyl groups (not in rings)
    # This could help exclude non-peptide amide-containing compounds
    amino = Chem.MolFromSmarts("[NH2,NH1;!$(NC=O)]")
    carboxyl = Chem.MolFromSmarts("[CX3](=O)[OH]")
    if not mol.HasSubstructMatch(amino) or not mol.HasSubstructMatch(carboxyl):
        return False, "Missing free amino or carboxyl groups"
    
    return True, "Contains multiple amide bonds and characteristics of a peptide antibiotic"
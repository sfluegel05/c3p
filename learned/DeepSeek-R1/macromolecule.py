"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a high molecular weight molecule with repeating structural units.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    
    # Adjusted molecular weight threshold for biopolymers
    if mol_wt < 500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low"

    # Protein detection: multiple amide bonds in backbone pattern
    # Looking for consecutive amide bonds (at least 10 for small proteins)
    amide_pattern = Chem.MolFromSmarts("[CX3][CX3](=[OX1])[NX3H]")
    amide_count = len(mol.GetSubstructMatches(amide_pattern))
    if amide_count >= 10:
        return True, f"Contains {amide_count} amide bonds, indicative of a protein"

    # Polysaccharide detection: multiple glycosidic linkages
    # Improved pattern for any oxygen connecting two carbons in rings (excluding esters)
    glycosidic_pattern = Chem.MolFromSmarts("[C&r][OX2;!$(O=C)][C&r]")
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))
    if glycosidic_count >= 6:
        return True, f"Contains {glycosidic_count} glycosidic-like linkages"

    # Check for repeating units using fingerprint complexity
    # Macromolecules tend to have high numbers of unique substructures
    unique_branched_substructs = rdMolDescriptors.CalcNumUniqueBranchedSubstructs(mol)
    if unique_branched_substructs > 15 and mol_wt > 1000:
        return True, "Complex structure with repeating units"

    # Fallback for very large molecules
    if mol_wt > 3000:
        return True, f"Very high molecular weight ({mol_wt:.1f} Da)"

    return False, "Does not meet macromolecule criteria"